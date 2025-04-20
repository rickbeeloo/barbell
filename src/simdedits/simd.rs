#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
use std::arch::x86_64::*;

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
pub const MAX_QUERIES: usize = 32;


#[cfg(target_arch = "aarch64")]
use std::arch::aarch64::*;

#[cfg(target_arch = "aarch64")]
pub const MAX_QUERIES: usize = 16;

// Define SIMD types per architecture
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
pub type SimdVector = __m256i;

#[cfg(target_arch = "aarch64")]
pub type SimdVector = uint8x16_t;

#[derive(Debug)]
pub struct TransposedQueries {
    pub vectors: Vec<SimdVector>,
    pub query_length: usize,
    pub n_queries: usize,
}

impl TransposedQueries {
    pub fn new(queries: Vec<&[u8]>) -> Self {
        assert!(!queries.is_empty(), "No queries provided");
        assert!(
            queries.iter().all(|q| q.len() == queries[0].len()),
            "All queries must have the same length"
        );

        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        {
            if is_x86_feature_detected!("avx2") {
                unsafe { Self::new_x86_avx2(queries) }
            } else {
                panic!("AVX2 support is required for x86/x86_64 architecture");
            }
        }

        #[cfg(target_arch = "aarch64")]
        {
            unsafe {
                return Self::new_aarch64_neon(queries);
            }
        }

        #[cfg(not(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64")))]
        {
            panic!("Unsupported architecture");
        }
    }

    // x86/x86_64 with AVX2
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    #[target_feature(enable = "avx2")]
    unsafe fn new_x86_avx2(queries: Vec<&[u8]>) -> Self {
        let query_length = queries[0].len();
        let n_queries = queries.len();
        assert!(n_queries <= 32, "Number of queries exceeds 32");

        let mut vectors = Vec::with_capacity(query_length);
        for i in 0..query_length {
            let mut query_bases = [0u8; 32];
            for (q, query) in queries.iter().enumerate() {
                if q < n_queries {
                    query_bases[q] = query[i];
                }
            }
            let v = _mm256_loadu_si256(query_bases.as_ptr() as *const __m256i);
            vectors.push(v);
        }

        Self {
            vectors,
            query_length,
            n_queries,
        }
    }

    // aarch64 with NEON
    #[cfg(target_arch = "aarch64")]
    #[target_feature(enable = "neon")]
    unsafe fn new_aarch64_neon(queries: Vec<&[u8]>) -> Self {
        let query_length = queries[0].len();
        let n_queries = queries.len();
        assert!(n_queries <= 16, "Number of queries exceeds 16");

        let mut vectors = Vec::with_capacity(query_length);
        for i in 0..query_length {
            let mut query_bases = [0u8; 16];
            for (q, query) in queries.iter().enumerate() {
                if q < n_queries {
                    query_bases[q] = query[i];
                }
            }
            let v = vld1q_u8(query_bases.as_ptr());
            vectors.push(v);
        }

        Self {
            vectors,
            query_length,
            n_queries,
        }
    }
}


pub fn simd_search(
    transposed: &TransposedQueries,
    target: &[u8],
) -> Vec<u8> {
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if is_x86_feature_detected!("avx2") {
            unsafe { align_x86_avx2(transposed, target) }
        } else {
            panic!("AVX2 support is required for x86/x86_64 architecture");
        }
    }
    
    #[cfg(target_arch = "aarch64")]
    {
        // On Apple Silicon (macOS), NEON is always available
        // so we can directly call the function
        #[cfg(target_os = "macos")]
        unsafe {
            return align_aarch64_neon(transposed, target);
        }
        
        // For other aarch64 platforms, check for NEON support
        #[cfg(not(target_os = "macos"))]
        {
            if std::arch::is_aarch64_feature_detected!("neon") {
                unsafe { align_aarch64_neon(transposed, target) }
            } else {
                panic!("NEON support is required for aarch64 architecture");
            }
        }
    }
    
    #[cfg(not(any(target_arch = "x86", target_arch = "x86_64", target_arch = "aarch64")))]
    {
        panic!("Unsupported architecture");
    }
}

// x86/x86_64 with AVX2
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn align_x86_avx2(
    transposed: &TransposedQueries,
    target: &[u8],
) -> Vec<u8> {
    let TransposedQueries { vectors: query_transposed, query_length, n_queries } = transposed;

    let target_length = target.len();
    let one = _mm256_set1_epi8(1);
    let row_size = target_length + 1;

    let mut dp = vec![_mm256_setzero_si256(); row_size * 2];
    let base_ptr = dp.as_mut_ptr();
    let mut min_vec = _mm256_set1_epi8(u8::MAX as i8);

    for i in 1..=*query_length {
        let band = (i & 1) as isize;
        let curr = base_ptr.offset(band * row_size as isize);
        let prev = base_ptr.offset((1 - band) * row_size as isize);
        curr.write(_mm256_set1_epi8(i as i8));

        let qv = query_transposed[i - 1];
        let mut j = 1;
        let unroll_bound = target_length.saturating_sub(3);

        while j <= unroll_bound {
            for k in 0..4 { 
                let tj = target[j - 1 + k];
                let tv = _mm256_set1_epi8(tj as i8);
                let eq = _mm256_cmpeq_epi8(qv, tv);
                let sub = _mm256_andnot_si256(eq, one);

                let d = _mm256_add_epi8(prev.add(j - 1 + k).read(), sub);
                let u = _mm256_add_epi8(prev.add(j + k).read(), one);
                let l = _mm256_add_epi8(curr.add(j - 1 + k).read(), one);
                let m1 = _mm256_min_epu8(d, u);
                let mc = _mm256_min_epu8(m1, l);

                curr.add(j + k).write(mc);
                if i == *query_length {
                    min_vec = _mm256_min_epu8(min_vec, mc);
                }
            }
            j += 4;
        }

        // Handle remaining elements
        while j <= target_length {
            let tv = _mm256_set1_epi8(target[j - 1] as i8);
            let eq = _mm256_cmpeq_epi8(qv, tv);
            let sub = _mm256_andnot_si256(eq, one);
            let d = _mm256_add_epi8(prev.add(j - 1).read(), sub);
            let u = _mm256_add_epi8(prev.add(j).read(), one);
            let l = _mm256_add_epi8(curr.add(j - 1).read(), one);
            let m1 = _mm256_min_epu8(d, u);
            let mc = _mm256_min_epu8(m1, l);
            curr.add(j).write(mc);
            if i == *query_length {
                min_vec = _mm256_min_epu8(min_vec, mc);
            }
            j += 1;
        }
    }

    let mut tmp = [0u8; 32];
    _mm256_storeu_si256(tmp.as_mut_ptr() as *mut __m256i, min_vec);
    tmp[..*n_queries].to_vec()
}

// aarch64 implementation with NEON
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn align_aarch64_neon(
    transposed: &TransposedQueries,
    target: &[u8],
) -> Vec<u8> {
    let TransposedQueries { vectors: query_transposed, query_length, n_queries } = transposed;

    let target_length = target.len();
    let one = vdupq_n_u8(1);
    let row_size = target_length + 1;

    let mut dp = vec![vdupq_n_u8(0); row_size * 2];
    let base_ptr = dp.as_mut_ptr();
    let mut min_vec = vdupq_n_u8(u8::MAX);

    for i in 1..=*query_length {
        let band = (i & 1) as isize;
        let curr = base_ptr.offset(band * row_size as isize);
        let prev = base_ptr.offset((1 - band) * row_size as isize);
        curr.write(vdupq_n_u8(i as u8));

        let qv = query_transposed[i - 1];
        let mut j = 1;
        let unroll_bound = target_length.saturating_sub(3);

        while j <= unroll_bound {
            for k in 0..4 {
                let tj = target[j - 1 + k];
                let tv = vdupq_n_u8(tj);
                let eq = vceqq_u8(qv, tv);
                // No AndNot like in avx2, so vmvnq and vandq
                let sub = vandq_u8(vmvnq_u8(eq), one);

                let d = vaddq_u8(prev.add(j - 1 + k).read(), sub);
                let u = vaddq_u8(prev.add(j + k).read(), one);
                let l = vaddq_u8(curr.add(j - 1 + k).read(), one);
                let m1 = vminq_u8(d, u);
                let mc = vminq_u8(m1, l);

                curr.add(j + k).write(mc);
                if i == *query_length {
                    min_vec = vminq_u8(min_vec, mc);
                }
            }
            j += 4;
        }

        // Handle remaining elements
        while j <= target_length {
            let tv = vdupq_n_u8(target[j - 1]);
            let eq = vceqq_u8(qv, tv);
            let sub = vandq_u8(vmvnq_u8(eq), one);
            let d = vaddq_u8(prev.add(j - 1).read(), sub);
            let u = vaddq_u8(prev.add(j).read(), one);
            let l = vaddq_u8(curr.add(j - 1).read(), one);
            let m1 = vminq_u8(d, u);
            let mc = vminq_u8(m1, l);
            curr.add(j).write(mc);
            if i == *query_length {
                min_vec = vminq_u8(min_vec, mc);
            }
            j += 1;
        }
    }

    let mut tmp = [0u8; 16];
    vst1q_u8(tmp.as_mut_ptr(), min_vec);
    tmp[..*n_queries].to_vec()
}

#[cfg(test)]
mod tests {
    use super::*;
    use pa_bitpacking::search::*;
    use std::time::Instant;
    use rand::Rng;

    #[test]
    fn test_all_zero() {
        let queries = vec![b"ACGT".as_slice(), b"ACGT".as_slice(), b"ACGT".as_slice()];
        let simd_query_set = TransposedQueries::new(queries);
        let target = b"ACGT";
        let edits = simd_search(&simd_query_set, target);
        assert_eq!(edits, vec![0, 0, 0]);
    }

    #[test]
    fn test_edits() {
        let queries = vec![b"ACGT".as_slice(), b"AGGT".as_slice(), b"ATTT".as_slice()];
        let simd_query_set = TransposedQueries::new(queries);
        let target = b"ACGT";
        let edits = simd_search(&simd_query_set, target);
        assert_eq!(edits, vec![0, 1, 2]);
    }

    #[test]
    fn more_complex_edits() {
        let queries = vec![b"AACTAGCTAGCATCGA".as_slice(), b"AGCTAGGTAGAATCGT".as_slice(), b"AGCAAAAAAGAATCGT".as_slice()];
        let simd_query_set = TransposedQueries::new(queries);
        let target = b"AGCTAGCTAGCATCGT";
        let edits = simd_search(&simd_query_set, target);
        assert_eq!(edits, vec![2, 2, 5]);
    }

    fn random_dna(len: usize) -> Vec<u8> {
        let mut rng = rand::thread_rng();
        let dna: Vec<u8> = (0..len).map(|_| {
            let c = ['A', 'T', 'C', 'G'][rng.gen_range(0..4)];
            c as u8
        }).collect();
        dna
    }

    
}