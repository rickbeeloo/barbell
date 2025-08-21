# Barbell - a new perspective on demultiplexing

### Why barbell?
- More accurate barcode scoring scheme, around **5x less 
barcode bleeding** than Dorado
- Over **100x less trimming errors** compared to Dorado
- Equivalent, or **better assemblies**
- **Contamination free** assemblie by removing artefact reads
- Easily applicable to **custom experiments**
- Still **very fast**


### The steps of a barbell workflow?
Barbell is a read demultiplexer and trimmer that focuses on detecting experimental errors. A general workflow follows four steps:
1. **Annotate**, which finds barcodes and their flanks in the reads. It uses [edit distance](https://www.biorxiv.org/content/10.1101/2025.07.22.666207v1) to locate the flanks and a new scoring scheme based on [Lodhi et al.](https://www.jmlr.org/papers/v2/lodhi02a.html) to score barcodes
2. **Inspect**, based on the annotate results Barbell can "print" all the reads patterns showing where barcodes are, if they are concatenated, etc. 
3. **Filter**, you specify what *patterns* are valid in your opinion and where to cut reads
4. **Trim** takes the filtered output and the original reads, trims them, and splits them based on the barcodes

## Installing Barbell
Barbell is written in Rust

#### Install Rust
Check if you have Rust installed, run `rustc --version`, if no output you should install Rust:

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

### Important note
**Note**: requires **AVX2 and BMI2** as these are used
in [sassy](https://github.com/RagnarGrootKoerkamp/sassy). It is therefore unlikely to work on Apple Silicon or ARM-based systems, and may also fail on older CPUs. <u>Most modern x86_64 Ubuntu systems should work, depending on the CPU</u>.

### From releases (when public)
see [releases](https://github.com/rickbeeloo/barbell/releases), for example:
```bash
wget https://github.com/rickbeeloo/barbell-sg/releases/download/v0.1.8/barbell-sassy-rewrite-x86_64-unknown-linux-gnu.tar.xz

```

### From source
```bash
git clone https://github.com/rickbeeloo/barbell-sg
cd barbell-sg
cargo build --release
```


The `barbell` executable is now in `/target/release` 

## Annotate
```bash
Usage: barbell annotate [OPTIONS] --input <INPUT> --queries <QUERIES>

Options:
  -i, --input <INPUT>
          Input FASTQ file
  -t, --threads <THREADS>
          Number of threads [default: 10]
  -o, --output <OUTPUT>
          Output file path [default: output.tsv]
  -q, --queries <QUERIES>
          Query files (comma-separated paths)
  -b, --barcode-types <BARCODE_TYPES>
          Barcode types (comma-separated: Fbar,Rbar,Fflank,Rflank) matching your query file (-q) [default: Fbar]
      --flank-max-errors <INT>
          Flank maximum erors in flank, ONLY set manually when you know what you are doing
      --verbose
          Enable verbose output for debugging
      --min-score <MIN_SCORE>
          Barcode: fraction compared to 'perfect' match score for top candidate [default: 0.5]
      --min-score-diff <MIN_SCORE_DIFF>
          Barcode: fraction difference between top 2 candidates [default: 0.05]
  -h, --help
          Print help
```

### Simple
The basic input for `annotate` is:
```bash 
barbell annotate -i reads.fastq -q queries.fasta
```

- `--input`, having your reads in fastq
- `--queries`, the sequences you are looking for (in fasta), this *should* have the flanks as well for example if you are demultiplexing a rapid barcoding kit you create a file like [this](examples/rapid_bars.fasta) where we have
`<left_flank><barcode><right_flank>` for Dorado these are mostly on their chemistry webpage of which most are also in [examples](examples/). 


In case you have mutliple queries, for example for dual end you can do:
```bash 
barbell -i reads.fastq -q query1.fasta,query2.fasta -b Fbar,Rbar
```
The `-b` options tells barbell that `query1.fasta` has forward barcodes, and `query2.fasta` reverse barcodes. This does not change the algorithm it just makes filtering easier.
For native barcoding the left and right sequences are the same so you could just run:
```bash
barbell -i reads.fastq -q query1.fasta -b Fbar
```
and just in the filtering step allow it to match in `rc` on the right side. 

### Extra options
It's often useful to control how "strict" barbell is find barcodes, this is determined by:
1. `--min-score <X>`, to be very safe you could use `0.6`, for assembly where it matters less around `0.4`, on average around `0.5`
2. `--min_score_diff <X>`, generally between `0.05-0.07` seems fine. 


### Advanced 
If your flanks are short for example in case of native barcoding the fitted edit distance formula will allow up to 2 edits in the flanks.
If you have too little reads you could try to up this using `--flank-max-errors 3` which overrides the automatic edit cut off.




### What are patterns?
At  the core of Barbell are its *patterns* that describe a read:
They always have the following format:
```
Tag[<orientation>, <label>, <rel_pos>, <cut>]
```
Where:
- Tag: `Ftag` or `Rtag`
- orientation: `fw` or `rc`, relative to your query input file
- label: 
     - `*` match any label from input fasta 
     - `?1`, wildcard, for example in dual barcoding you can use `Ftag[fw, ?1, @left(0,250)]__Ftag[rc, ?1, @right(0,250)]`. Can use as many wildcards as you want, `?2`, `?3`, etc.
     - `~g1`, substring, 


