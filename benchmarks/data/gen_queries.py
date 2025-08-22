import fastapy

rapid_left = "GCTTGGGTGTTTAACC"
rapid_right = "GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA"

# Read rapid bars as fasta
with open("rapid_queries.fasta", "w") as q1, open("../../examples/rapid_bars.fasta", "w") as q2:
    for record in fastapy.parse("rapid_bars.txt"):
        barcode_name = record.id
        barcode_seq = record.seq
        # Create query
        query = rapid_left + barcode_seq + rapid_right
        # Create record
        q1.write(f">{barcode_name}\n{query}\n")
        q2.write(f">{barcode_name}\n{query}\n")

native_left = "ATTGCTAAGGTTAA"
native_right = "CAGCACCT"
 
def reverse_complement(seq):
    # use make trans atgc 
    trans_table = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C"
    }
    return "".join(trans_table[base] for base in seq[::-1])

with open("native_queries.fasta", "w") as q1, open("../../examples/native_bars.fasta", "w") as q2:
    for record in fastapy.parse("native_bars.txt"):
        barcode_name = record.id
        barcode_seq = record.seq
        barcode_seq = reverse_complement(barcode_seq)
        # Create query
        query = native_left + barcode_seq + native_right
        # Create record
        q1.write(f">{barcode_name}\n{query}\n")
        q2.write(f">{barcode_name}\n{query}\n")

