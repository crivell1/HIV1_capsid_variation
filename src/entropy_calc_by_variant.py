import sys
import numpy as np
import pandas as pd

df = pd.read_csv(sys.argv[1])

# Calculate contributions only for valid, positive counts.
df["entropy_value"] = np.nan
valid = (df["count"] > 0) & (df["valid_sequences"] > 0)

p = df.loc[valid, "count"] / df.loc[valid, "valid_sequences"]
df.loc[valid, "entropy_value"] = -p * np.log2(p)

entropy = (
    df.groupby(["variant", "alignment_column"])
    .agg(
        entropy_bits=(
            "entropy_value",
            lambda values: values.sum(min_count=1)
        ),
        valid_sequences=("valid_sequences", "first"),
        group_sequences=("group_sequences", "first")
    )
    .reset_index()
)

entropy.to_csv(sys.argv[2], index=False)