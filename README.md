# Barbell Demultiplexer

## Overview
Barbell demultiplexer uses a [custom semi-global aligner](https://github.com/RagnarGrootKoerkamp/astar-pairwise-aligner) designed to accurately identify barcodes, even in read ends where semi-global constraints typically introduce artificial edits. Barbell operates with fewer parameters and is optimized to detect as much contamination as possible.

The process begins by generating edit distances for every position in the read. Minima are then selected from these distances and backtraced to obtain full alignments. All obtained alignments compete for the same read region, with the best alignment ultimately being selected.

Barbell's focus extends beyond barcode identification to include flanking sequences such as adapters. In cases where concatenated barcodes are nearly undetectable, these flanking sequences may still be present. Identifying such cases is crucial to prevent misannotations downstream.

## Features
📜 Annotate FASTQ files with barcode information
🚀 Supports multi-threading for faster processing
⚙️ Optional autotuning for selecting the main parameter
🧩 Modular design with a flexible strategy for demultiplexing


## Installation

To install Barbell Demultiplexer, clone the repository and build it using Cargo:

```sh
# Clone the repository
git clone git@github.com:rickbeeloo/barbell-sg.git
cd barbell-sg

# Build the project
RUSTFLAGS='-C target-cpu=native' cargo build --release
```

Note, we use compile flags (`RUSTFLAGS='-C target-cpu=native'`) to target the native architecture as the aligner uses [SIMD instructions](https://en.wikipedia.org/wiki/SIMD).

## Usage

### Annotate FASTQ Files
To annotate a FASTQ file with barcode information, use the `annotate` command:
Note, the executable is in `./target/release/`

```sh
./barbell annotate -i input.fastq -o output.txt -q queries.fasta -t 8 --tune
```

**Options:**
- `-i, --input` (required): Input FASTQ file
- `-o, --output` (required): Output file path
- `-q, --queries` (required): Query files (comma-separated paths)
- `-t, --threads` (optional): Number of threads (default: 5)
- `--tune` (optional): Enable autotuning

It's recommended to use the `--tune` flag to find the best parameter for your input queries unless you already have a good estimamte. This is simply the fraction of the input sequence which is allowed to be mutated. Generally two random sequences have an edit distance of 0.5, so this cut off should be below that, by default `0.35`.

You can provide multiple query files, so for example if you have a dual-end barcode experiment you can provide `left.fasta` and `right.fasta`, which will  be prefixed wiht `L` and `R` respectively in the output. An example query file, for rapid barcoding can be found in [examples/rapid_barcodes.fasta](examples/rapid_barcodes.fasta).


### Plot Results
(Not implemented yet)

```sh
./barbell plot
```

## Example Output
```sh
Starting annotation...
Processing input.fastq with 8 threads...
Annotation complete!
```

## License
This project is licensed under the MIT License.

## Contributing
Contributions are welcome! Feel free to submit a pull request or open an issue.

## Contact
For any questions or issues, please open an issue on GitHub

