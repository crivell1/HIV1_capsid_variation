import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("p24_entropy.csv")
plt.plot(df["alignment_column"], df["entropy_bits"])
plt.xlabel("Alignment Column")
plt.ylabel("Entropy (bits)")
plt.title("Entropy Across Alignment Columns")
plt.show()

