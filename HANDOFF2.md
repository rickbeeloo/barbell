## 1) Removed `--adapter-trim` from `barbell trim` and `barbell kit`

- Deleted the four adapter flags (`--adapter-trim`, `--adapter-min-match-len`, `--adapter-min-identity`, `--adapter-end-slack`) from both subcommands in `bin/main.rs`. They now fail loudly with `error: unexpected argument '--adapter-trim' found`.
- Removed `adapter_trim: AdapterTrimConfig` from `TrimConfig` and `KitConfig`, and the plumbing in `src/kits/use_kit.rs` and `src/trim/trim.rs`.
- Reverted `process_read_and_anno` to its pre-adapter signature and dropped the now-obsolete `test_adapter_trim_shifts_barcode_cuts` test and the `AdapterTrim`/`AdapterTrimmer` import from `trim.rs`.
- Removed the now-meaningless `enabled` field from `AdapterTrimConfig` (adapter trimming is now defined purely by the command).
- Docs updated: README section now explains the two-step `kit` then `trim-adapters` workflow, CHANGELOG and HANDOFF reflect the new state.

## 2) `barbell trim-adapters` is now multithreaded

- New `-t/--threads` option (default 10, `DEFAULT_ADAPTER_THREADS`).
- `trim_adapters_in_fastq` reads paraseq batches of 512 records, trims each batch via rayon `par_chunks_mut(...).for_each_init(...)` (each rayon worker gets its own `AdapterTrimmer` via `for_each_init`, so the adapter isn't re-encoded per read), then writes the batch **in input order** — output is deterministic and independent of thread count. Multiple input files are also processed concurrently, one file per task.
- `threads <= 1` uses a plain sequential loop (no thread pool built at all).
- Verified: 200k reads / 900 MB → 116 s (1 thread), 33 s (4), 18 s (12), with **identical md5** across all three, and gzip + multi-file runs likewise identical.
- Tests: 66 pass; new `test_trim_adapters_in_fastq_multithreaded_matches_single_threaded` (byte-identical 1 vs 8 threads) and an extended multi-file test comparing 1 vs 4 threads. `cargo clippy` shows no new warnings vs. baseline; `cargo fmt` clean for the touched files.
