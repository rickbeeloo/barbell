# 🦀 Barbell Demultiplexer (Rust)

## Overview
Barbell demultiplexer uses a [custom semi-global aligner](https://github.com/RagnarGrootKoerkamp/astar-pairwise-aligner) designed to accurately identify barcodes, even in read ends where semi-global constraints typically introduce artificial edits. Barbell operates with fewer parameters and is optimized to detect as much contamination as possible.

The process begins by generating edit distances for every position in the read. Minima are then selected from these distances and backtraced to obtain full alignments. All obtained alignments compete for the same read region, with the best alignment ultimately being selected.

Barbell's focus extends beyond barcode identification to include flanking sequences such as adapters. In cases where concatenated barcodes are nearly undetectable, these flanking sequences may still be present. Identifying such cases is crucial to prevent misannotations downstream.

## Installation (CLI)

To install Barbell Demultiplexer, clone the repository and build it using Cargo:

```sh
# Clone the repository
git clone git@github.com:rickbeeloo/barbell-sg.git
cd barbell-sg

# Build the project
RUSTFLAGS='-C target-cpu=native' cargo build --release
```

If you don't have Rust installed, you can install it using (see [Rust installation guide](https://www.rust-lang.org/tools/install)):

```sh
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

Note, we use compile flags (`RUSTFLAGS='-C target-cpu=native'`) to target the native architecture as the aligner uses [SIMD instructions](https://en.wikipedia.org/wiki/SIMD).

## Quickstart
Inside the folder `target/release/`:
```
./barbell annotate -i input.fastq -o output.txt -q queries.fasta -t 8 --tune
```

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

You can provide multiple query files, so for example if you have a dual-end barcode experiment you can provide `left.fasta` and `right.fasta`, which will  be prefixed wiht `L` and `R` respectively in the output. _DONT_ combine your left and right queries in the same file as this will fail to automatically extract the shared flanking regions. An example query file, for rapid barcoding can be found in [examples/rapid_barcodes.fasta](examples/rapid_barcodes.fasta).

## Features
- 📜 Annotate FASTQ files with barcode information
- 🚀 Supports multi-threading for faster processing
- ⚙️ Optional autotuning for selecting the main parameter
- 🧩 Modular design with a flexible strategy for demultiplexing


### Plot Results
(Not implemented yet)

```sh
./barbell plot
```

## Example Output
The output file will have the following format:
```sh
read	label	start	end	edits	dist.to.end	read.len
acc308e2-e1b8-49f6-bc33-ba38b4b1facc	1R#F_fw	6	124	5	6	3627
b3e54c72-c72d-4553-bee9-b07245a723f8	11R#F_fw	0	115	10	1	3452
b3e54c72-c72d-4553-bee9-b07245a723f8	1F#F_fw	3351	3452	7	-101	3452
110aa41d-b36e-46cc-af12-51bade6b7a2b	1F#R_rc	1	113	4	1	3807
110aa41d-b36e-46cc-af12-51bade6b7a2b	1R#R_rc	3712	3807	21	-95	3807
ea7a7b9a-fa35-493f-8b0f-5e3dba93ba32	1F#R_rc	5	122	2	5	3424
ea7a7b9a-fa35-493f-8b0f-5e3dba93ba32	5R#R_rc	3326	3424	12	-98	3424
eeef1710-5c0f-4720-85f1-821fd4cca1b0	2R#F_fw	6	118	15	6	3395
eeef1710-5c0f-4720-85f1-821fd4cca1b0	1F#F_fw	3295	3395	11	-100	3395
1f38871f-88d1-442e-a0a8-0aa1305ab4cc	3R#F_fw	9	122	3	9	3432
1f38871f-88d1-442e-a0a8-0aa1305ab4cc	1F#F_fw	3351	3432	22	-81	3432
3fc6772e-3638-4fef-b764-b230aa1feec0	1F#R_rc	8	120	7	8	3417
3fc6772e-3638-4fef-b764-b230aa1feec0	7R#R_rc	3343	3417	25	-74	3417
a98c36bd-5b77-49ab-af05-f6d113b8ba7e	1F#R_rc	8	123	3	8	3377
a98c36bd-5b77-49ab-af05-f6d113b8ba7e	2R#R_rc	3280	3377	9	-97	3377
a9d1a073-3af3-45ed-9984-35e5b58764d7	1F#R_rc	3	117	2	3	3430
a9d1a073-3af3-45ed-9984-35e5b58764d7	1F#F_fw	3354	3430	30	-76	3430
bc6da060-b185-4723-aeb8-5d90015321fa	4R#F_fw	0	115	6	1	3943
bc6da060-b185-4723-aeb8-5d90015321fa	1F#F_fw	3866	3943	20	-77	3943
```

With the read identifier (`read`), followed by the `label` - the most important column. It will have the sequence identifiers from your input file with as suffix `#F_rw`, `#F_rc`, `#R_fw` and `#R_rc` for forward, forward rc, reverse and reverse rc matches respectively. It's always smart to take a look at the file to understand what was matched, and in what orientation and whether this is what you expect. The label column can also contain `Flank` this indicates that no barcode was good enough to be aligned but that the flank was still clearly present. The `edits` is the number of edits compared to the query. The `dist.to.end` is the distance to the end of the read which is positive when closer to the left end of the read and negative when closer to the right end. This makes it easy to filter later on for example saying your barcodes should be at least 250bp away from the read ends. The read length column is the length of the read.

## License
This project is licensed under the MIT License.

## Contributing
Contributions are welcome! Feel free to submit a pull request or open an issue.

## Contact
For any questions or issues, please open an issue on GitHub

