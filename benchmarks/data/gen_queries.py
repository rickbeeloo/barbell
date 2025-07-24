import fastapy

left = "GCTTGGGTGTTTAACC"
right = "GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA"

# Read rapid bars as fasta
with open("rapid_queries.fasta", "w") as q:
    for record in fastapy.parse("rapid_bars.txt"):
        barcode_name = record.id
        barcode_seq = record.seq
        # Create query
        query = left + barcode_seq + right
        # Create record
        q.write(f">{barcode_name}\n{query}\n")





