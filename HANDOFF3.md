Done. Here's what I found.

## No — the batch size isn't the problem

I instrumented the pipeline (`trim_adapter_file` phase timing + per-thread `/proc` sampling) and measured it directly:

| Check                                        | Result                                               |
|----------------------------------------------|------------------------------------------------------|
| Time in `trim` (parallel) vs fill/copy/write | **98.6%** vs <1.5%                                   |
| Per-worker activity (10 threads)             | every worker running, evenly balanced (~10% each)    |
| Pure compute scaling (no I/O)                | 1t=5.8k → 12t=42.5k reads/s = **7.3× at 12 threads** |
| **Uncompressed** input, `-t 10`              | **844% CPU**                                         |
So `ADAPTER_BATCH_RECORDS = 512` is not the limiter — work is distributed evenly and compute scales almost linearly. I left it alone.

## What actually caused the low CPU

Two serial, single-core stages sit outside the rayon pool:

1. **Gzip *input* decode** — `paraseq` → `niffler` → `flate2`/`miniz_oxide`, entirely on one thread. It pins a core while the trim threads wait, capping utilization (it was ~10% of wall on realistically-compressible data, ~2% on highly compressible FASTQ).
2. **Gzip *output* encode** — the `flate2::write::GzEncoder` is also single-threaded. This is the dominant cost when writing `.gz`: 44s vs 13s for the same data, with CPU dropping to ~290%.

## What I did: rapidgzip-rust for input decode

`rapidgzip-core 0.3.1` exposes a `Read + Send` reader that plugs straight into paraseq (it's literally built for that). New `open_batch_reader()`:
- **gzip regular files** → rapidgzip parallel decode, using the `-t/--threads` budget
- **plain regular files** → plain `BufReader` (1 MiB) — *faster* than before, since the old path paid a `niffler` sniff + flate2 type-erasure
- **FIFOs / `/dev/stdin`** → unchanged niffler path (parallel paths need positional reads)

I found and fixed a real bug while doing this: my first version sniffed the gzip magic bytes on the *same* handle, consuming 2 bytes and corrupting plain input. It now probes a throwaway handle, with extension fallback.

**Results** (500 MB, 34k long ONT-like reads):

|                               | before          | after                 |
|-------------------------------|-----------------|-----------------------|
| gzip in → plain out, `-t 10`  | 14.3s, 699% CPU | **13.4s, 756% CPU**   |
| plain in → plain out, `-t 10` | 12.7s           | 12.8s (no regression) |
| gzip in → gzip out, `-t 10`   | 44.3s           | 43.0s, ~290% CPU      |
Output is **byte-identical** (md5-matched) across 1 vs 10 threads, gzip vs plain input, plain vs gzip output, FIFO input, gzip content stored under a `.fastq` name, and `/dev/stdin`. 68 tests pass (2 new), `cargo fmt`/`clippy` clean, README + CHANGELOG updated.

## Two caveats

- **The gain here is modest** (14.3s → 13.4s) because on FASTQ, deflate decode is genuinely cheap. The payoff is larger on less-compressible input or slow storage, and it removes the "one core is pinned while 9 idle" failure mode. Decode is now ~7× faster in isolation (0.27s vs 1.98s on a 500 MB payload).
- **Gzip *output* is still single-threaded** and is now the clear bottleneck for `.gz` output — that's a larger change needing a parallel deflate crate. If you want, I can look at that next; it's what would actually lift CPU utilization on your `--gzip` runs.
