import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv(sys.argv[1])

## had to look up how to create a heatmap with matplotlib - apparently its not as straightforward as using another heatmap library like seaborn
## a bit of googling and HopGPT here to get this working correctly

heatmap = df.pivot( # have to pivot the dataframe to create a matrix suitable for a heatmap
    index="variant",
    columns="alignment_column",
    values="entropy_bits"
).sort_index(axis=1)

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

ticks = np.arange(0, len(heatmap.columns), 20)
ax.set_xticks(ticks)
ax.set_xticklabels(heatmap.columns[ticks])

ax.set_xlabel("Alignment column")
ax.set_ylabel("Variant")
ax.set_title("Shannon entropy by HIV-1 variant")
fig.colorbar(image, ax=ax, label="Entropy (bits)")

fig.tight_layout()
fig.savefig(sys.argv[2], dpi=150)
plt.show()