import pandas as pd
import numpy as np

input_csv = sys.argv[1] #continue with data/p24_entropy_by_variant.csv

entropy_df = pd.read_csv(input_csv)

# Create a pivot table with variants as rows and alignment columns as columns
pivot_df = entropy_df.pivot(index="variant", columns="alignment_column", values="entropy_bits")

annotation_df = pd.read_csv("data/HX2B_annotations_by_aa.tsv", sep="\t")
annotation_df = annotation_df.rename(columns={"amino_acid": "alignment_column"})

pivot_df = pivot_df.merge(annotation_df, left_index=True, right_on="alignment_column", how="left")

## this gives us a dataframe with the alignment column, the entropy bits for each variant, and the annotation for that alignment column

