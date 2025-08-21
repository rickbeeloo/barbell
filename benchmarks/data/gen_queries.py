import fastapy

rapid_left = "GCTTGGGTGTTTAACC"
rapid_right = "GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA"

# Read rapid bars as fasta
with open("rapid_queries.fasta", "w") as q:
    for record in fastapy.parse("rapid_bars.txt"):
        barcode_name = record.id
        barcode_seq = record.seq
        # Create query
        query = rapid_left + barcode_seq + rapid_right
        # Create record
        q.write(f">{barcode_name}\n{query}\n")

native_left = "ATTGCTAAGGTTAA"
native_right = "CAGCACCT"
 
with open("native_queries.fasta", "w") as q:
    for record in fastapy.parse("native_bars.txt"):
        barcode_name = record.id
        barcode_seq = record.seq
        # Create query
        query = native_left + barcode_seq + native_right
        # Create record
        q.write(f">{barcode_name}\n{query}\n")

