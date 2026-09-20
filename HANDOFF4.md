## What changed

`barbell trim-adapters --gzip` previously wrote output through a single-threaded `flate2::write::GzEncoder`. HANDOFF3 already identified this as the remaining bottleneck ("Gzip *output* is still single-threaded… that's a larger change needing a parallel deflate crate"), which is exactly what I addressed.

The fix, in `src/trim/adapter.rs`:

- Added a `gzip_writer()` helper that wraps the output file in [`gzp`](https://github.com/sstadick/gzp)'s `ZBuilder::<Gzip>`, deflating blocks on a background thread pool sized by `-t/--threads`.
- Added a `gzip_threads` field to `AdapterJob`, wired from `config.threads`, so each output writer gets the requested budget.
- `threads == 1` falls through to a plain synchronous writer (gzp's `ZBuilder` does this automatically), so single-threaded runs are unaffected.
- Added `gzp` with the pure-Rust `deflate_rust` backend (no C/zlib-ng dependency), matching the existing pure-Rust `flate2`/miniz_oxide setup.

This is the write-side counterpart to the existing parallel `rapidgzip-core` read path, so `-t/--threads` now covers reading, trimming, **and** writing.

## Verification

- **71 tests pass** (3 new): multithreaded gzip output matches single-threaded; `.gz` suffix vs `--gzip` flag behavior; parallel-compress → parallel-decompress round trip.
- **No new clippy/fmt warnings** in the touched file.
- `cargo fmt` had reformatted two unrelated pre-existing files, so I reverted those to keep the diff focused.

End-to-end on a 114 MB / 40k-read ONT-like FASTQ:

| Scenario                     | Before             | After              |
|------------------------------|--------------------|--------------------|
| Plain in → gzip out, `-t 10` | 20.4s, **99% CPU** | **3.2s, 946% CPU** |
| gzip in → gzip out, `-t 10`  | 8.9s, **350% CPU** | **3.2s, 947% CPU** |
Output is a valid single-member gzip file (`gzip -t` passes) that decompresses to **byte-identical** reads (md5-matched) across 1 vs 10 threads, plain vs gzip input, and single vs multi-file folder mode. Empty input and the `.gz`-content-under-a-`.fastq`-name case both work.

One note: multi-file folder mode names its outputs `*.adapter_trimmed.fastq` (no `.gz` suffix even when gzip is on) — that's pre-existing naming behavior in `plan_jobs`, unrelated to this change, so I left it alone. Let me know if you'd like the suffix added too.
