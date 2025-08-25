# Barbell - a new perspective on demultiplexing

### Why barbell?
- More accurate barcode scoring scheme, around **2-5x less 
barcode bleeding** than Dorado
- Over **100x less trimming errors** compared to Dorado
- Equivalent, or **better assemblies**
- **Contamination free** assembly by removing artefact reads
- Easily applicable to **custom experiments**
- Still **very fast**

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


## Quickstart
Most Nanopore kits (except those for RNA and Twist) were added to the `kit` command and can be used without 
much manual interference: (if you have a custom kit [jump here[(#todo)]). 
While these "presets" help you to quickly analyze data and get a feel of barbell we highly recommend 
going through the "Understanding Barbell" section as you will likely get much more from your reads 
by understanding them. 

```bash
barbell kit -k <kit-name> -i <reads.fastq> -o <output_folder> --maximize
```
The `--maximize` option is generally recommended (e.g. for assembly) unless you really want the near perfect barcode configuration. 

### Native barcoding example (SQK-NBD114-96)
```bash
barbell kit -k SQK-NBD114-96 -i <reads.fastq> -o <output_folder> --maximize
```
This use our lower bound edit distance formula to set a cut-off for the sequences flanking the barocde. 
This is generally quite strict but works well for good sequence quality. If you notice a lot of reads are missed during the annotate step you could try to lower it to for example `--flank-max-errors 5`, but make sure you again check the inspect output to see if you are not creating too many "Fflank" matches which likely indicate random matches, thereby 
decreasing your reads.

### Rapid barcoding example (SQK-RBK114-96)
```bash
barbell kit -k SQK-RBK114-96 -i <reads.fastq> -o <output_folder> --maximize
```
The currently supported kits can be found [here[(data/supported_kits.txt)]









