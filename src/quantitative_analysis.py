import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
import sys
from pathlib import Path

aa_entropy_df = pd.read_csv(sys.argv[1])
aa_entropy_variant_df = pd.read_csv(sys.argv[2])
output_dir = Path(sys.argv[3])
output_dir.mkdir(parents=True, exist_ok=True) ## pathlib stuff

## need to analyze both csvs to see whether there's any correlation between strucutral annotation and entropy, with covariates of SASA and variant

required_columns = {"aa_position", "entropy_bits", "major_region", "minor_region", "rsa"}

def summarize_regions(df, region_column, by_variant = False):
    group_columns = [region_column]
    if by_variant == True:
        group_columns = ["variant", region_column]
    summary = df.groupby(group_columns).agg(
        positions = ("entropy_bits", "count"),
        mean_entropy=("entropy_bits","mean"),
        median_entropy=("entropy_bits","median"),
        std_entropy=("entropy_bits","std")
    )
    return summary.reset_index() # had to add the .reset_index() to 


# run functions and get results into the results/ directory

major_summary = summarize_regions(aa_entropy_df, "major_region")
minor_summary = summarize_regions(aa_entropy_df, "minor_region")

variant_major_summary = summarize_regions(
    aa_entropy_variant_df, "major_region", by_variant=True
)

variant_minor_summary = summarize_regions(
    aa_entropy_variant_df, "minor_region", by_variant=True
)

# Add regression function here - CAR
# 
# 
# 
# 



major_summary.to_csv(output_dir / "pooled_major_regions.csv", index=False)
minor_summary.to_csv(output_dir / "pooled_minor_regions.csv", index=False)
variant_major_summary.to_csv(output_dir / "variant_major_regions.csv", index=False)
variant_minor_summary.to_csv(output_dir / "variant_minor_regions.csv", index=False)

pooled_correlation.to_frame().T.to_csv(
    output_dir / "pooled_rsa_entropy_correlation.csv", index=False
)
pooled_variant_correlations.to_csv(
    output_dir / "variant_rsa_entropy_correlations.csv", index=False
)
print(f"Saved six analysis tables to {output_dir}")
