Everything works. Here's a summary of what I implemented.

## What was added

An option to trim Oxford Nanopore **kit14 ligation adapters** (the AMX114 adapters used by `SQK-LSK114`, `SQK-ULK114`, and all other "114" kits) from read ends, in the spirit of Porechop / `dorado trim`.

### Why this needed a new mechanism
Barbell's existing trimming is *barcode-driven*: reads are cut at the barcode. But in ligation libraries the adapter sits **outside** the barcode (`[adapter][barcode cassette][insert]`), so it survives barcode trimming. There was no "adapter-only" search in the codebase (this is exactly the long-standing request in issue #5, where users specifically mention LSK114/ULK114/AMX114).

### Adapter sequences (authoritative, from Dorado)
- Front adapter: `CCTGTACTTCGTTCAGTTACGTATTGC`
- Reverse complement (other capture strand): `GCAATACGTAACTGAACGAAGTACAGG`

### New code
**`src/trim/adapter.rs`** — `AdapterTrimmer`, built on the existing `sassy` matcher. It:
- searches both adapter orientations with partial-alignment support (basecalling often truncates the first bases of a read),
- only accepts adapters **anchored at the very read end** (middle adapters are deliberately left alone — those need splitting, not end trimming),
- **chains tandem adapters** (a real ligation artefact, e.g. `[adapter][adapter][insert]`),
- tolerates a few noisy bases in front of the adapter,
- ships a standalone `trim_adapters_in_fastq` that preserves read order/orientation and handles plain/gzipped FASTQ.

### Integration
- New standalone `barbell trim-adapters` command. Adapter trimming is *only*
  available here (the `--adapter-trim` flags on `barbell trim` / `barbell kit`
  were removed again).
- `barbell trim-adapters` is multithreaded (`-t/--threads`, default 10) and
  preserves read order within each file even with several threads.
- New `AdapterTrimConfig`, plus tuning flags `--min-match-len`, `--min-identity`, `--end-slack`.

### Validation
- **66 tests pass** (22 new), including a real LSK114/LSK110 read from Dorado's test suite, tandem adapters, partial adapters, middle-adapter rejection, and that 8-thread output is byte-identical to 1-thread output.
- Tuned on 6,000 simulated kit14 reads + 12,000 negatives: **99.0% sensitivity, 0.13% false-positive rate** on random sequence and near-zero on reads with only internal adapters.
- Verified end-to-end: a simulated `SQK-RBK114-24` dataset demultiplexed into correct barcodes with 0/240 reads retaining adapter or barcode sequence, and reverse-complemented reads handled correctly.
- Measured on a 200k-read/900MB FASTQ: 116s with 1 thread, 33s with 4, 18s with 12 — all three byte-identical.
- `cargo build --release`, `cargo test`, and `cargo clippy` are clean (no new warnings vs. baseline).

Docs were added to the README (with a Quick-link entry) and a CHANGELOG note referencing issue #5.
