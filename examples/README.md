

# 🦀 Barbell Demultiplexer  - Examples 

## Rapid barcode demultiplexing

### 1. Annotate

We query our reads, using competitive alignments, for rapid barcodes and their flanks (see this file for the rapid barcodes: [rapid_barcodes.fasta](rapid_barcodes.fasta)).
```sh
./barbell annotate -i reads.fastq -o annotations.txt -q rapid_barcodes.fasta
```

### 2. Filter

We now filter for valid rapid barcode matches, the default filter you probably want is:
```sh
./barbell filter -i annotations.txt -o filtered.txt -p "Fbarcode[fw, *, >>,  @left(0 to 250)]"
```
We did notice that sometimes the same barcode is present twice at the left side, since we can still d determine  the sample we can include this in the filter (see [rapid_filters.txt](rapid_filters.txt)). You can also use the file using the `-f` options instead o f the `-p` option:

```sh
./barbell filter -i annotations.txt -o filtered.txt -f rapid_filters.txt
```

In the `rapid_filters.txt` we use `@left(-250 to 250)` as semi-global alignments can overlap. 

### 3. Trim

We now trim the reads based on the filtered annotation file. In rapid barcoding we only really care about what barcode they had as we just have a `fw` barcode, so orientation is not important. We can exclude orientation from the output file names using `--no-orientation`.


```sh
./barbell trim -i filtered.txt -r reads.fastq -o trimmed --no-orientation
```

That's it!  


## Native barcode demultiplexing
todo!




