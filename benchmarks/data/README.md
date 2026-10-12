### Simulated data evals
We generate randoms sequences for the following categories:
- GroupI: Just random data - no flanks nor barcodes are present
- GroupII: Perfect example, we insert a flank + barcode on the left (as prefix)
- GroupIII: Same as groupII except that now we randomly trim up to 10nts from left or right of the sequecne
- GroupIV: GroupII + extra barcode + flank pair directly after first with 0..10 gap
- GroupV: GroupII + barcode + flank in the middle of the read 
- GroupVI: GroupII + reverse_complement(flank + barcode) at the right end of the read

So for this we expect to not find any matches for `GroupI`, all matches for `GroupII`, as many as possible 
for `GroupIII` and "nothing" for `GroupIV-VI` as these should be considered invalid reads for the 
rapid barcoding protocol.


### Rapid barcode data
The used barcodes are not clearly mentioned on the website but are in the code:

```C++
// RBK_1_96 is the same as BC_1_96 except for 26, 39, 40, 58, 54 and 60.
const std::vector<std::string> RBK_1_96 = {
        "BC01", "BC02", "BC03", "BC04",  "BC05",  "BC06",  "BC07",  "BC08", "BC09", "BC10",  "BC11",
        "BC12", "BC13", "BC14", "BC15",  "BC16",  "BC17",  "BC18",  "BC19", "BC20", "BC21",  "BC22",
        "BC23", "BC24", "BC25", "RBK26", "BC27",  "BC28",  "BC29",  "BC30", "BC31", "BC32",  "BC33",
        "BC34", "BC35", "BC36", "BC37",  "BC38",  "RBK39", "RBK40", "BC41", "BC42", "BC43",  "BC44",
        "BC45", "BC46", "BC47", "RBK48", "BC49",  "BC50",  "BC51",  "BC52", "BC53", "RBK54", "BC55",
        "BC56", "BC57", "BC58", "BC59",  "RBK60", "BC61",  "BC62",  "BC63", "BC64", "BC65",  "BC66",
        "BC67", "BC68", "BC69", "BC70",  "BC71",  "BC72",  "BC73",  "BC74", "BC75", "BC76",  "BC77",
        "BC78", "BC79", "BC80", "BC81",  "BC82",  "BC83",  "BC84",  "BC85", "BC86", "BC87",  "BC88",
        "BC89", "BC90", "BC91", "BC92",  "BC93",  "BC94",  "BC95",  "BC96"};
```

Though with a mistake, should be "48" different instead of "58". For example "RBK26" is `ACTATGCCTTTCCGTGAAACAGTT`
which is not mentioned anywhere on the [chemistry doc](https://nanoporetech.com/document/chemistry-technical-document#barcode-sequences). So we compile the list `rapid_bars.txt` based on their code.

