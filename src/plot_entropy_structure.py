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
    chain_id = "H" # FreeSASA notation for chain A -> the entire protein is in chain A

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

def plot_entropy_sasa_scatter(merged_df, output_file):
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

def line_plot_entropy_sasa(merged_df, output_file):

    plot_df = merged_df.sort_values("aa_position").copy()
    fig, axes = plt.subplots(
        2, 1, figsize=(14, 7), sharex=True
    )
    axes[0].plot(plot_df["aa_position"], plot_df["entropy_bits"], color="black")
    axes[1].plot(plot_df["aa_position"], plot_df["rsa"], color="steelblue")
    axes[0].set_ylabel("Entropy (bits)")
    axes[1].set_ylabel("Relative solvent accessibility")
    axes[1].set_xlabel("HXB2 p24 amino-acid position")




    ## STOP here
    plt.figure(figsize=(10, 6))
    sns.lineplot(data=merged_df, x="aa_position", y="entropy_bits") #, hue="variant")
    
    # Set plot labels and title
    plt.xlabel("Solvent Accessible Surface Area (SASA)")
    plt.ylabel("Shannon Entropy (bits)")
    plt.title("Entropy vs SASA") # by Variant")
    
    # Save the plot to a file
    plt.savefig(output_file)
    plt.close()
    
if __name__ == "__main__":
    
    #print(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])

    annotated_entropy_df = load_annotated_entropy(sys.argv[1], sys.argv[2])
    print(annotated_entropy_df.head())  # Debugging line to check the structure of the annotated entropy DataFrame

    sasa_df = calculate_sasa(sys.argv[3])
    print(sasa_df.head())  # Debugging line to check the structure of the SASA DataFrame

    merged_df = merge_sasa_entropy(annotated_entropy_df, sasa_df)
    print(merged_df.head())  # Debugging line to check the structure of the merged DataFrame

    plot_entropy_sasa_scatter(merged_df, sys.argv[4])
    line_plot_entropy_sasa(merged_df, sys.argv[5])

    # Lot of system arguments but trying to keep it modular and flexible


