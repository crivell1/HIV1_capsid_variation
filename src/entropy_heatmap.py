import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO

df = pd.read_csv(sys.argv[1])

# taking this from entropy_viz.py and modifying it to create a heatmap of entropy by variant instead of a line plot of entropy by amino acid position
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
df = df.merge(
    mapping_df,
    on="alignment_column",
    how="left",
    validate="many_to_one"
)
df = df.dropna(subset=["aa_position"])
df["aa_position"] = df["aa_position"].astype(int)

## had to look up how to create a heatmap with matplotlib - apparently its not as straightforward as using another heatmap library like seaborn
## a bit of googling and HopGPT here to get this working correctly

letter_variants = df["variant"].str.fullmatch(
    r"(?:A[1-8]?|B|C|D|F[12]?|G|H|J|K|L|M|N|O|P|U)+",
    na=False
)
df = df.loc[letter_variants].copy()
minimum_sequences = 1
df = df.loc[df["group_sequences"] >= minimum_sequences].copy()
if df.empty:
    raise ValueError("No variants meet the label and sequence-count filters.")

heatmap = df.pivot( # have to pivot the dataframe to create a matrix suitable for a heatmap
    index="variant",
    columns="aa_position",
    values="entropy_bits"
).sort_index(axis=1)

subtype_order = [
    "A", "B", "C", "D", "F", "G", "H",
    "J", "K", "L", "M", "N", "O", "P", "U"
]
order_lookup = {subtype: position for position, subtype in enumerate(subtype_order)}
variant_order = sorted(
    heatmap.index,
    key=lambda variant: (order_lookup.get(variant[0], 999), len(variant), variant)
)
heatmap = heatmap.reindex(variant_order)


group_sizes = df.groupby("variant")["group_sequences"].first() #group each variant and get the first group size for labeling
labels = [
    f"{variant} (n={int(group_sizes[variant])})"
    for variant in heatmap.index
]

cmap = plt.get_cmap("viridis").copy() # create a copy of the colormap so I can modify the missing vals without affecting the original
cmap.set_bad("lightgray")

fig, ax = plt.subplots(
    figsize=(16, max(3, len(heatmap.index) * 0.5))
)

image = ax.imshow(
    np.ma.masked_invalid(heatmap.to_numpy()), # hiding the missing values in the heatmap
    aspect="auto",
    interpolation="nearest",
    cmap=cmap,
    vmin=0,
    vmax=np.log2(20)
)

ax.set_yticks(range(len(labels))) ## axis formatting and labeling - will have to continue messing with this to get the same graph as the paper we are recreating
ax.set_yticklabels(labels)

row_variants = heatmap.index.tolist()
for row_number in range(1, len(row_variants)):
    previous_letter = row_variants[row_number - 1][0]
    current_letter = row_variants[row_number][0]
    if current_letter != previous_letter:
        ax.axhline(
            y=row_number - 0.5,
            color="white",
            linewidth=2
        )

ticks = np.arange(0, len(heatmap.columns), 20)
ax.set_xticks(ticks)
ax.set_xticklabels(heatmap.columns[ticks])

ax.set_xlabel("HXB2 p24 Amino Acid Position")
ax.set_ylabel("Variant")
ax.set_title("Shannon entropy by HIV-1 variant")
fig.colorbar(image, ax=ax, label="Entropy (bits)")

fig.tight_layout()
fig.savefig(sys.argv[3], dpi=150)
plt.show() #uncomment this line if you want to display the plot interactively
plt.close()