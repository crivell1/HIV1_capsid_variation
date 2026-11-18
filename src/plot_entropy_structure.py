import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import freesasa

def load_annotated_entropy(entropy_file, annotation_file):
    # Load the entropy and annotation data
    entropy_df = pd.read_csv(entropy_file)
    annotation_df = pd.read_csv(annotation_file, sep="\t")
    
    # Merge the two DataFrames on amino acid position
    merged_df = pd.merge(entropy_df, 
                         annotation_df, 
                         left_on="aa_position", 
                         right_on="amino_acid", 
                         how="left", 
                         validate="one_to_one") 
    
    # Drop the redundant amino acid column from the annotation DataFrame
    merged_df = merged_df.drop(columns=["amino_acid"])
    
    return merged_df
## this gives us a dataframe with the alignment column, the entropy bits for each variant, and the annotation for that alignment column

def calculate_sasa(pdb_file):
    # Load the SASA data
    structure = freesasa.Structure(pdb_file)
    result = freesasa.calc(structure)
    chain_id = "A" # FreeSASA notation for chain A -> the entire protein is in chain A

    # print(result.residueAreas())  # Debugging line to check the structure of the result
    # print(result.residueAreas.keys())  # Debugging line to check the available chain IDs
    # print(result.residueAreas()[chain_id])  # Debugging line to check the structure of the result

    sasa_records = []

    for res_num_str, area in result.residueAreas()[chain_id].items():
    # Strip whitespace and keep only standard integer residue numbers
        clean_num = res_num_str.strip()
        if clean_num.isdigit():
            sasa_records.append({
                "aa_position": int(clean_num),
                "sasa": area.total,
                "rsa": area.relativeTotal if area.hasRelativeAreas else None
            })

    sasa_df = pd.DataFrame(sasa_records)

    return sasa_df

def merge_sasa_entropy(annotated_entropy_df, sasa_df):

    # Merge the annotated entropy DataFrame with the SASA DataFrame on amino acid position
    merged_df = pd.merge(annotated_entropy_df, 
                         sasa_df, 
                         left_on="aa_position", 
                         right_on="aa_position", 
                         how="left", 
                         validate="one_to_one")
    
    # Drop the redundant amino acid column from the SASA DataFrame
    #merged_df = merged_df.drop(columns=["amino_acid"])
    
    return merged_df

def plot_entropy_sasa_scatter_major(merged_df, output_file):
    # Create a scatter plot of entropy vs SASA
    plt.figure(figsize=(10, 6))
    sns.scatterplot(data=merged_df, x="rsa", y="entropy_bits", hue="major_region")
    
    # Set plot labels and title
    plt.xlabel("Solvent Accessible Surface Area (SASA)")
    plt.ylabel("Shannon Entropy (bits)")
    plt.title("Entropy vs SASA") # by Variant")
    
    # Save the plot to a file
    plt.savefig(output_file)
    plt.close()

def plot_entropy_sasa_scatter_minor(merged_df, output_file):
    # Create a scatter plot of entropy vs SASA
    plt.figure(figsize=(10, 6))
    sns.scatterplot(data=merged_df, x="rsa", y="entropy_bits", hue="minor_region")
    
    # Set plot labels and title
    plt.xlabel("Solvent Accessible Surface Area (SASA)")
    plt.ylabel("Shannon Entropy (bits)")
    plt.title("Entropy vs SASA") # by Variant")
    
    # Save the plot to a file
    plt.savefig(output_file)
    plt.close()

def line_plot_entropy_sasa_major(merged_df, output_file):

    #Had to look up all of this matplotlib code to get the backgroun shading for regions. 

    plot_df = merged_df.sort_values("aa_position").copy()
    fig, axes = plt.subplots(
        2, 1, figsize=(14, 7), sharex=True
    )

    # Creating two plots, one for entropy and one for SASA, with the same x-axis (amino acid position)
    axes[0].plot(plot_df["aa_position"], plot_df["entropy_bits"], color="black")
    axes[1].plot(plot_df["aa_position"], plot_df["rsa"], color="steelblue")

    #labeling the plots and axes
    axes[0].set_ylabel("Entropy (bits)")
    axes[1].set_ylabel("Relative solvent accessibility")
    axes[1].set_xlabel("HXB2 p24 amino-acid position")

    # shading of each strucutral region  
    region_colors = {"NTD": "#cfe8ff", "IDR": "#ffe6ad", "CTD": "#d9eed3"}
    region_changes = plot_df["major_region"].ne(plot_df["major_region"].shift())
    region_blocks = region_changes.cumsum()

    # region shading and labeling for each structural region
    for _, block in plot_df.groupby(region_blocks): # _ is a placeholder to iterate through the groupby object without needing the group name
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

    # Save it and move on to the next plot
    fig.tight_layout()
    fig.savefig(output_file, dpi=150)
    plt.close(fig)

def line_plot_entropy_sasa_minor(merged_df, output_file):

# almost identical to the major region plot, but with minor regions instead of major regions.
    plot_df = merged_df.sort_values("aa_position").copy()
    fig, axes = plt.subplots(
        2, 1, figsize=(14, 7), sharex=True
    )
    axes[0].plot(plot_df["aa_position"], plot_df["entropy_bits"], color="black")
    axes[1].plot(plot_df["aa_position"], plot_df["rsa"], color="steelblue")
    axes[0].set_ylabel("Entropy (bits)")
    axes[1].set_ylabel("Relative solvent accessibility")
    axes[1].set_xlabel("HXB2 p24 amino-acid position")

    # Keys match the minor_region labels in HXB2_annotations_by_aa.tsv 
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
        
    fig.tight_layout()
    fig.savefig(output_file, dpi=150)
    plt.close(fig)
    
if __name__ == "__main__":
    
    #print(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])

    annotated_entropy_df = load_annotated_entropy(sys.argv[1], sys.argv[2])
    # print(annotated_entropy_df.head())  # Debugging line to check the structure of the annotated entropy DataFrame

    sasa_df = calculate_sasa(sys.argv[3])
    # print(sasa_df.head())  # Debugging line to check the structure of the SASA DataFrame

    merged_df = merge_sasa_entropy(annotated_entropy_df, sasa_df)
    # print(merged_df.head())  # Debugging line to check the structure of the merged DataFrame
    merged_df.to_csv("data/p24_entropy_aa_with_annotations.csv")

    plot_entropy_sasa_scatter_major(merged_df, sys.argv[4])
    plot_entropy_sasa_scatter_minor(merged_df, sys.argv[5])
    line_plot_entropy_sasa_major(merged_df, sys.argv[6])
    line_plot_entropy_sasa_minor(merged_df, sys.argv[7])

    # Lot of system arguments but trying to keep it modular and flexible


