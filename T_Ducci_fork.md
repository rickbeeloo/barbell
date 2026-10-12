# T_Ducci fork of barbell

## Goal

Reattach methylation (MM/ML) to reads after trimming, outside of barbell.
This requires the cut coordinates on the original read, which upstream barbell does not output.

## `--return-cuts-idxs` flag

Available in both `barbell trim` and `barbell kit`. Off by default, in which case the output is identical to upstream barbell.

When enabled, three TAB-separated SAM tags are appended to the header of each output fragment:

```
@read1 MM:Z:...	ML:B:C,...	ch:i:12	bs:i:60	be:i:900	bf:i:0
```

| Tag | Meaning |
|---|---|
| `bs:i:` | fragment start on the original read (0-based) |
| `be:i:` | fragment end on the original read (exclusive) |
| `bf:i:` | `1` if the fragment was reverse complemented (`--flip`), `0` otherwise. Always present |

- **Coordinates**: always refer to the original read and are computed before flipping. The fragment is `orig_seq[bs:be]`, reverse complemented if `bf = 1`.
- **`--skip-trim`**: the whole read is written, but `bs`/`be` still report where the cut would have happened.
- **Chimeric reads**: each fragment (`read1`, `read1_1`, …) gets its own coordinates.
- **Header**: barbell's format (space between read id and original tags) is unchanged. Checked with samtools 1.21 that `samtools import -T '*'` correctly parses both the original tags and `bs`/`be`/`bf`.

## Changes

- `bin/main.rs`: flag added to `Commands::Trim` and `Commands::Kit`, propagated to `TrimConfig` / `KitConfig`.
- `src/config.rs`: `return_cuts_idxs` field added to `TrimConfig` and `KitConfig`.
- `src/kits/use_kit.rs`: flag passed from `KitConfig` to `TrimConfig`.
- `src/trim/trim.rs`:
  - `process_read_and_anno` takes the new `return_cuts_idxs` parameter and returns a fifth element, `Option<(start, end, flipped)>`, which is `None` when the flag is off.
  - `trim_matches` passes the flag through and writes the tags to the header.
  - Tests: updated existing calls and added assertions on the coordinates for the plain, `skip_trim` and `flip` cases.
