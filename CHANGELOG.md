### Changes in new release

<!-- next-header -->

## git

## 0.3.2

## 0.3.2

- Added:
  - Gzip support thanks @Vanishborn
  - Add option to skip trimming (#17)
- Fixed:
  - Prevent progress overwriting errors, such as failing to open enough file handles (`ulimit`)
  - Update cigar parsing to new sassy version (2c42eb3)

## Earlier

- Use new version of Sassy with Mac support (v1.4)
- Set default cut off to `--min-score=0.2` and `--min-score-diff=0.1` — loose but good for assembly
- Use "kits" using `barbell kit --<kit-name>`