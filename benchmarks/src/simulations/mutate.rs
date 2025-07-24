use rand::Rng;
use std::collections::HashMap;

/// Mutate sequence with at most max_edits - from sassy repo
pub fn mutate_sequence(sequence: &[u8], min_edits: usize, max_edits: usize) -> Vec<u8> {
    let mut rng = rand::rng();
    let mut seq = sequence.to_vec();
    for _ in 0..rng.random_range(min_edits..=max_edits) {
        let idx = rng.random_range(0..seq.len());
        match rng.random_range(0..3) {
            0 => {
                let current = seq[idx];
                let mut new_char;
                // Keep trying until we get a different character
                loop {
                    new_char = b"ACGT"[rng.random_range(0..4)];
                    if new_char != current {
                        break;
                    }
                }
                seq[idx] = new_char;
            }
            1 if seq.len() > 1 => {
                seq.remove(idx);
            }
            2 => seq.insert(idx, b"ACGT"[rng.random_range(0..4)]),
            _ => {}
        }
    }
    seq
}

pub fn random_trim_side(seq: &[u8], max_trim: usize) -> Vec<u8> {
    let mut rng = rand::rng();
    let mut seq = String::from_utf8_lossy(seq).to_string();
    let trim_length = rng.random_range(1..=max_trim);
    let trim_front = rng.random_bool(0.5);
    let trim_back = rng.random_bool(0.5);

    if trim_front {
        seq.drain(0..trim_length);
    }

    if trim_back {
        seq.drain(seq.len() - trim_length..seq.len());
    }

    seq.as_bytes().to_vec()
}

#[cfg(test)]
mod tests {

    #[test]
    fn test_random_trim_side() {
        let seq: &'static [u8; 12] = b"ACGTACGTACGT";
        let trimmed = super::random_trim_side(seq, 2);
        assert!(trimmed.len() < seq.len());
    }

    #[test]
    fn test_mutate_sequence() {
        let seq: &'static [u8; 12] = b"ACGTACGTACGT";
        let mutated = super::mutate_sequence(seq, 1, 2);
        assert!(mutated != seq);
    }
}
