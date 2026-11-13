import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO

df = pd.read_csv(sys.argv[1])
#alignment colum included the new blank regions from the altered p24 pull

alignment_file = sys.argv[2]  # Read the mapping file

hxb2_seq = None

for record in SeqIO.parse(alignment_file, "fasta"):
    if record.id.split(".")[-1] == "K03455":
        hxb2_seq = str(record.seq)
if hxb2_seq is None:
    raise ValueError("HXB2 reference K03455 was not found.")

mapping_rows = []
aa_position =0 
for alignment_column, amino_acid in enumerate(hxb2_seq, start =1): #had to get some google help on the enumerate function to get the alignment column and amino acid
    if amino_acid != "-" and amino_acid != ".":
        aa_position += 1
        mapping_rows.append({"alignment_column": alignment_column, "aa_position": aa_position})
    else:
        mapping_rows.append({"alignment_column": alignment_column, "aa_position": None})


mapping_df = pd.DataFrame(mapping_rows, columns=["alignment_column", "aa_position"])
# Merge the mapping DataFrame with the entropy DataFrame on alignment_column
df = df.merge(mapping_df, on="alignment_column", how="left")
df = df.dropna(subset=["aa_position"])  # Drop rows where aa_position is NaN
df.sort_values(by="aa_position", inplace=True)  # Sort by aa_position
plt.plot(df["aa_position"], df["entropy_bits"])
plt.xlabel("Amino Acid Position")
plt.ylabel("Entropy (bits)")
plt.title("Entropy Across Amino Acid Position of HXB2 p24")
plt.show() # uncomment this line if you want to display the plot interactively
plt.savefig(sys.argv[3], dpi=150)
plt.close()

