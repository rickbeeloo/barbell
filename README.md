# 🦀 Barbell Demultiplexer (Rust)

## Overview
Barbell demultiplexer uses a [custom semi-global aligner](https://github.com/RagnarGrootKoerkamp/astar-pairwise-aligner) designed to accurately identify barcodes, even in read ends where semi-global constraints typically introduce artificial edits. Barbell operates with fewer parameters and is optimized to detect as much contamination as possible. Moreover, Barbell is not restricted to any specific experimental setup and can detect a variety of complex barcode patterns. Barbell follows three steps:

1. <b> Annotate </b>, where are your queries? What is their score? Is the barcode clear?
2. <b> Filter </b>, what reads match the experimental design? `Fbar--Rbar`, or `Fbar--`, etc.
3. <b> Trim </b>, Cut of all contamination, barcodes, adapters, etc.

Enjoy!

## Installation (CLI)

### Easy - from releases
The easiest way to install Barbell is from the [releases page](https://github.com/rickbeeloo/barbell-sg/releases).
1. Download
2. Unzip
3. use executable🚀

<i>Note since these binaries are meant to  be stable across architectures the aligner is slower than building from source.</i>


### From source

To install Barbell Demultiplexer, clone the repository and build it using Cargo. You  have to swtich to the nigthly channel if you are not already:

```sh
# Clone the repository
git clone git@github.com:rickbeeloo/barbell-sg.git
cd barbell-sg

# Can skip if already on nightly
rustup install nightly
rustup override set nightly

# Build the project
RUSTFLAGS='-C target-cpu=native' cargo build --release
```

If you don't have Rust installed, you can install it using (see [Rust installation guide](https://www.rust-lang.org/tools/install)):

```sh
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

Note, we use compile flags (`RUSTFLAGS='-C target-cpu=native'`) to target the native architecture as the aligner uses [SIMD instructions](https://en.wikipedia.org/wiki/SIMD). <b>This will likely give better performance than using the binaries from the releases page.</b>

## Quickstart
Inside the folder `target/release/`:
```
./barbell annotate -i input.fastq -o output.txt -q queries.fasta -t 8 --tune
```

## Usage

Barbell follows a  three step process:
1. Annotate, find all querie sequences in the reads 
2. Filter, filter for relevant patterns (e.g. fbar-----rbar)
3. Trim, use results from `Filter` to trim the reads 

### Annotate: FASTQ Files
To annotate a FASTQ file with barcode information, use the `annotate` command:
Note, the executable is in `./target/release/`

```sh
./barbell annotate -i reads.fastq -o annotations.txt -q queries.fasta -t 8 --tune
```

**Options:**
- `-i, --input` (required): Input FASTQ file
- `-o, --output` (required): Output file path
- `-q, --queries` (required): Query files (comma-separated paths)
- `-t, --threads` (optional): Number of threads (default: 5)
- `--tune` (optional): Enable autotuning

It's recommended to use the `--tune` flag to find the best parameter for your input queries unless you already have a good estimamte. This is simply the fraction of the input sequence which is allowed to be mutated. Generally two random sequences have an edit distance of 0.5, so this cut off should be below that, by default `0.35`.

You can provide multiple query files, so for example if you have a dual-end barcode experiment you can provide `left.fasta` and `right.fasta`, which will  be prefixed wiht `L` and `R` respectively in the output
- Provide your sequences looking at the top strand (not really a problem as reverse complement is also checked)
-  _DONT_ combine your left and right queries in the same file as this will fail to automatically extract the shared flanking regions.
An example query file, for rapid barcoding can be found in [examples/rapid_barcodes.fasta](examples/rapid_barcodes.fasta). 


### Filter: read annotation file

```sh
./barbell filter -i annotations.txt -o filtered.txt -p "Fbarcode[rc, *, >>, @left(0 to 250)]"
```

**Options:**
- `-i, --input` (required): Input file path
- `-o, --output` (required): Output file path
- `-p, --pattern` (required): Pattern string to filter by
- `-f, --file` (optional): File containing patterns to filter by

For quick checks or simple patterns, you can directly pass the pattern as a string using the `-p` option. For more complex filtering or when using multiple patterns, you can place them in a file and use the `-f` option instead. See [examples/rapid_filters.txt](examples/rapid_filters.txt) for example patterns. The patterns can be as complex as you want, see below for mmore information on building patterns.


### Trim: trim reads based on pattern results

<b> Note, we now keep the longest end which might not always be correct if the read is very short. This will be fixed in the next update</b>

```sh
./barbell trimm -i filtered.txt -r reads.fastq -o trimmed.fastq
```

**Options:**
- `-i, --input` (required): Input filtered annotation file
- `-r, --reads` (required): Input FASTQ file containing the reads
- `-o, --output` (required): Output FASTQ file path for trimmed reads
- `--no-label`: Disable label in output filenames
- `--no-orientation`: Disable orientation in output filenames
- `--no-flanks`: Disable flank in output filenames

This will trim reads based on the trim position markers in the pattern (`>>` or `<<`). By default, the output filenames will include labels, orientations, and flanks. You can disable these components in the output filenames using the corresponding flags. For exammple without any disable flags the files in the `output` folder will look like:
```bash
1R_fw__1F_fw.trimmed.fastq
2R_fw__1F_fw.trimmed.fastq
```

With orientation disabled:
```bash
1F__1R.trimmed.fastq
1F__2R.trimmed.fastq
```

In many default experiments you will mostly care about the presence of the barcode, and not the orientation, however it's recommended to look at the orientation information here (or in the filtered output file) to see if the orientation is correct.

## Understanding Patterns

Patterns are used to identify and filter specific sequences within reads. They help in specifying the expected structure of sequences, such as the presence of barcodes, their orientation, and their position within the read.

### Components of a Pattern

1. **Label**: Represents the identifier from the fastq file. It can be specific (e.g., `10R`) or a wildcard (`*` for any label, `?1` for matching labels).

2. **Orientation**: Indicates the direction of the sequence.
   - `fw`: Forward orientation (same as in the query file).
   - `rc`: Reverse complement orientation.

3. **Group**: Specifies the group the query belongs to.
   - `Fbarcode`: Forward barcode (from the first query file).
   - `Rbarcode`: Reverse barcode (from the second query file).
   - `Flank`: Flanking sequence (e.g., primer or adapter with an unknown barcode).

4. **Position**: Defines the location of the sequence within the read.
   - `@left(x to y)`: Distance from the left end of the read.
   - `@right(x to y)`: Distance from the right end of the read.
   - `@prev_left(x to y)`: Relative position to the previous match.

5. **Trimming**: Specifies the trimming position.
   - `>>`: Cut after the match.
   - `<<`: Cut before the match.

### Writing Patterns

#### Basic Pattern Structure

A pattern is written in the format: `Fbarcode[<extra_filters>]`, and multiple pattern elements are combined using `__`, for exammple `Fbarcode[<extra_filters>]__Rbarcode[<extra_filters>]`.

#### Filtering Barcodes

- **Left Side of the Read**:
  - Pattern: `Fbarcode[fw, *, @left(0 to 250)]`
  - Meaning: Look for `Fbarcode` in forward orientation, any label, within 0 to 250 bases from the left end.

- **Right Side of the Read**:
  - Pattern: `Rbarcode[rc, *, @right(0 to 250)]`
  - Meaning: Look for `Rbarcode` in reverse complement orientation, any label, within 0 to 250 bases from the right end.

#### Relative Positioning

- **Relative to Previous Match**:
  - Pattern: `Fbarcode[fw, ?1, @left(0 to 250)]__Fbarcode[fw, ?1, @prev_left(-100 to 200)]`
  - Meaning: The first `Fbarcode` is within 0 to 250 bases from the left end. The second `Fbarcode` is within -100 to 200 bases from the first match, and both must have the same label (`?1`).

#### Wildcards and Exact Matches

- **Specific and Any Barcode**:
  ```bash
  Rbarcode[fw, *, @left(0 to 250)]__Fbarcode[fw, 1F, @right(0 to 250)]
  Fbarcode[rc, 1F, @left(0 to 250)]__Rbarcode[rc, *, @right(0 to 250)]
  ```
  - Meaning: `1F` is a specific barcode, while `*` allows any barcode. The patterns account for both forward and reverse complement orientations.

- **Matching Labels**:
  ```bash
  Fbarcode[rc, ?1, @left(0 to 250)]__Rbarcode[rc, ?1, @right(0 to 250)]
  ```
  - Meaning: Both barcodes must have the same label (`?1`), useful for dual barcode experiments.

#### Trimming

- **Trimming Patterns**:
  ```bash
  Fbarcode[rc, ?1, >>, @left(0 to 250)]__Rbarcode[rc, ?1, <<, @right(0 to 250)]
  ```
  - Meaning: Cut after the `Fbarcode` match and before the `Rbarcode` match. `>>` trims at the end of the previous match, and `<<` trims at the start of the next match.

### Summary

- **Labels**: Use specific labels or wildcards (`*`, `?1`, etc.).
- **Orientation**: Specify `fw` or `rc`.
- **Groups**: Use `Fbar`, `Rbar`, or `Flank`.
- **Position**: Define using `@left`, `@right`, or `@prev_left`.
- **Trimming**: Use `>>` or `<<` for cutting sequences.

By understanding these components, you can create patterns tailored to your specific experimental needs. If you have any further questions or need additional clarification, feel free to ask!



## Features
- 📜 Annotate FASTQ files with barcode information
- 🚀 Supports multi-threading for faster processing 
- ⚙️ Optional autotuning for selecting the main parameter
- 🧩 Extremely modular and advanced filtering for custom experiments can be done using patterns

## License
This project is licensed under the MIT License.

## Contributing
Contributions are welcome! Feel free to submit a pull request or open an issue.

## Contact
For any questions or issues, please open an issue on GitHub

