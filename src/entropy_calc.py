import numpy as np
import pandas as pd
import sys

df =pd.read_csv(sys.argv[1])

p = df["count"] / df["valid_sequences"]
df["entropy_value"] = -p * np.log2(p)

entropy = (df.groupby("alignment_column")["entropy_value"].sum()
           .rename("entropy_bits").reset_index()) ## had to add this because the alignment positions didn't save, had to google this
entropy.to_csv(sys.argv[2], index=False)
