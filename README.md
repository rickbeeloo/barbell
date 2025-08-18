## Barbell - a new perspective on demultiplexing

### What is barbell?
Barbell is a read demultiplexer and trimmer that focuses on detecting experimental errors. A general workflow follows four steps:
1. **Annotate**, which finds barcodes and their flanks in the reads. It uses [edit distance](https://www.biorxiv.org/content/10.1101/2025.07.22.666207v1) to locate the flanks and a new scoring scheme based on [Lodhi et al.](https://www.jmlr.org/papers/v2/lodhi02a.html) to score barcodes
2. **Inspect**, based on the annotate results Barbell can "print" all the reads patterns showing where barcodes are, if they are concatenated, etc. 
3. **Filter**, you specify what *patterns* are valid in your opinion and where to cut reads
4. **Trim** takes the filtered output and the original reads, trims them, and splits them based on the barcodes

### Installing Barbell
Barbell is written in Rust

#### Install Rust
Check if you have Rust installed, run `rustc --version`, if no output you should install Rust:

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

#### Clone and build
```bash
git clone https://github.com/rickbeeloo/barbell-sg
cd barbell-sg
cargo build --release
```
The `barbell` executable is now in `/target/release` 


### What are patterns?
At  the core of Barbell are its *patterns* that describe a read:

