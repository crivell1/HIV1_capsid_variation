from plot_entropy_structure import merge_sasa_entropy, calculate_sasa, load_annotated_entropy
import matplotlib.pyplot as plt
import sys
import numpy as np
import pandas as pd
from Bio import SeqIO

### 
# stealing this all from the plot_entropy_structure.py script to get the annotated entropy dataframe and then merge it with the sasa dataframe to get a final dataframe that can be used for plotting
###

df = pd.read_csv(sys.argv[1])
#alignment colum included the new blank regions from the altered p24 pull

alignment_file = sys.argv[2]  # Read the mapping file "best_alignment.fasta"

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

# Had to add this to save the merged DataFrame to a new CSV file so that it can be used in the merge annotations script downstream
df["aa_position"] = df["aa_position"].astype(int)
df.to_csv(sys.argv[3], index=False)  # Save the merged DataFrame to a new CSV file

###
# Stealing from entropy_heatmap.py to filter the dataframe to only include letter variants and a minimum number of sequences
###

letter_variants = df["variant"].str.fullmatch(
    r"(?:A[1-8]?|B|C|D|F[12]?|G|H|J|K|L|M|N|O|P|U)+",
    na=False
)
df = df.loc[letter_variants].copy()
minimum_sequences = 1
df = df.loc[df["group_sequences"] >= minimum_sequences].copy()
if df.empty:
    raise ValueError("No variants meet the filters.")
annotation_df = pd.read_csv(sys.argv[4], sep="\t")

# Merge the annotated dataframe with the entropy dataframe
annotated_df = df.merge(
    annotation_df,
    left_on="aa_position", right_on="amino_acid",
    how="left", validate="many_to_one"
).drop(columns=["amino_acid"])

# calulate the sasa dataframe from the PDB file and merge it with the annotated dataframe
sasa_df = calculate_sasa(sys.argv[5])

merged_df = annotated_df.merge(
    sasa_df,
    on="aa_position",
    how="left", validate="many_to_one")

merged_df.to_csv("data/p24_variant_entropy_aa_with_annotations.csv")

###
# Stealing this function from plot_entropy_structure.py to create a line plot of entropy and sasa - Changing it to include the by variant lines
###

def line_plot_entropy_sasa_variant(merged_df, output_file, major_minor="major"):

    #Had to look up all of this matplotlib code to get the backgroun shading for regions. 

    plot_df = merged_df.sort_values("aa_position").copy()
    fig, axes = plt.subplots(
        2, 1, figsize=(14, 7), sharex=True
    )
    # sort alphabetically by variant to get a consistent color mapping across plots
    variants = sorted(plot_df["variant"].unique())
    colors = plt.get_cmap("turbo", len(variants))

    # for each variant, plot the entropy on the two subplots
    for index, variant in enumerate(variants):
        variant_df = plot_df.loc[plot_df["variant"] == variant]
        axes[0].plot(
            variant_df["aa_position"], variant_df["entropy_bits"],
            color=colors(index), label=variant, linewidth=0.8
        )
    # plot the sasa on the second subplot
    reference_df = plot_df.drop_duplicates(subset=["aa_position"])
    axes[1].plot(reference_df["aa_position"], reference_df["rsa"], color="steelblue")

    #labeling the plots and axes
    axes[0].set_ylabel("Entropy (bits)")
    axes[1].set_ylabel("Relative solvent accessibility")
    axes[1].set_xlabel("HXB2 p24 amino-acid position")

    # Keeping the region shading from the original plot_entropy_structure.py script - this is a bit of a hack but it works for now
    if major_minor == "major":
        region_colors = {"NTD": "#cfe8ff", "IDR": "#ffe6ad", "CTD": "#d9eed3"}
        region_changes = reference_df["major_region"].ne(reference_df["major_region"].shift())
        region_blocks = region_changes.cumsum()

        #Same as entropy_structure.py
        for _, block in reference_df.groupby(region_blocks): # _ is a placeholder to iterate through the groupby object without needing the group name
            region = block["major_region"].iloc[0]
            start = block["aa_position"].min() - 0.5
            end = block["aa_position"].max() + 0.5
            axes[0].text(
                            (start + end) / 2,
                            1.02,
                            region,
                            transform=axes[0].get_xaxis_transform(),
                            ha="center",
                            va="bottom",
                            fontsize=8,
                            rotation=90
                        )
            for ax in axes:
                ax.axvspan(start, end, color=region_colors.get(region, "lightgray"),
                        alpha=0.35, zorder=0)

    if major_minor == "minor":
        region_colors = {
                "B-hairpin": "#8dd3c7",
                "loop": "#eeeeee",
                "H1": "#a6cee3",
                "H2": "#b2df8a",
                "H3": "#fdbf6f",
                "H4": "#cab2d6",
                "H5": "#fb9a99",
                "H6": "#ffff99",
                "H7": "#80b1d3",
                "H8": "#b3de69",
                "H9": "#fccde5",
                "H10": "#bebada",
                "H11": "#ccebc5",
                "IDR": "#ffcc80",
                "MHR": "#d9b38c",
            }
        region_changes = plot_df["minor_region"].ne(plot_df["minor_region"].shift())
        region_blocks = region_changes.cumsum()
        
        for _, block in plot_df.groupby(region_blocks): # _ is a placeholder to iterate through the groupby object without needing the group name
            region = block["minor_region"].iloc[0]
            start = block["aa_position"].min() - 0.5
            end = block["aa_position"].max() + 0.5
            if region != "loop": # Only label non-loop regions to avoid cluttering the plot
                axes[0].text(
                            (start + end) / 2,
                            1.02,
                            region,
                            transform=axes[0].get_xaxis_transform(),
                            ha="center",
                            va="bottom",
                            fontsize=8,                        
                            rotation=90
                            )
            for ax in axes:
                ax.axvspan(start, end, color=region_colors.get(region, "lightgray"),
                            alpha=0.35, zorder=0)
    # Adding a legend for the variants to the first subplot
    axes[0].legend(
        title="Variant",
        loc="upper left",
        bbox_to_anchor=(1.02, 1),
        fontsize=7,
        ncol=2
    )

    #Save it!
    fig.tight_layout()
    fig.savefig(output_file, dpi=150, bbox_inches="tight")
    plt.close(fig)

if __name__ == "__main__":
    # sys arguments will be:
        # 1: variant entropy CSV
        # 2: final alignment FASTA
        # 3: mapped variant CSV output
        # 4: annotation TSV
        # 5: PDB file
        # 6: plot output

    # the rest of the code isn't written into functions since I was just copying and pasting from my previous scripts
    # call me lazy - that's fine with me

    line_plot_entropy_sasa_variant(merged_df, sys.argv[6], "major")
    line_plot_entropy_sasa_variant(merged_df, sys.argv[7], "minor")