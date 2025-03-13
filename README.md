# Barbell Demultiplexer

## Overview
Barbell demultiplexer uses a [custom semi-global aligner](https://github.com/RagnarGrootKoerkamp/astar-pairwise-aligner) that allows accurate barcode identification even at read ends where semi-global constrains normmally introduce artifical edits. In addition, Barbell uses less parameters and focusses on finding as much contamation as possible. It generates edit distances for every position in the read, from which minima are selected and backtraced to obtain full alignments... more explanation needed.

## Features
- Annotate FASTQ files with barcode information
- Supports multi-threading for faster processing
- Optional autotuning for performance optimization
- Modular design with a flexible strategy for demultiplexing

## Installation

To install Barbell Demultiplexer, clone the repository and build it using Cargo:

```sh
# Clone the repository
git clone git@github.com:rickbeeloo/barbell-sg.git
cd barbell-sg

# Build the project
RUSTFLAGS='-C target-cpu=native' cargo build --release
```

Note, we use compile flags to target the native architecture as the aligner uses [SIMD instructions](https://en.wikipedia.org/wiki/SIMD).

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

It's recommended to use the `--tune` flag to find the best parameter for your input queries unless you already have a good estimamte. This is simply the fraction of the input sequence which is allowed to be mutated. Generally two random sequences have an edit distnace of 0.5, so this cut off should be below that, by default `0.35`.

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

