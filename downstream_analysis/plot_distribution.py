import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("all_p_insertions.tsv", sep="\t", header=None)
# column indices depend on your TE_INFOS.bed format
# typically: chrom, start, end, TE_name, RS, freq, ...

freq_col = df.iloc[:, 10]   # adjust index to your actual freq column
rs_col   = df.iloc[:, 5]   # adjust index

fig, axes = plt.subplots(1,2, figsize=(10,4))
axes[0].hist(freq_col, bins=50)
axes[0].set_xlabel("Insertion frequency (%)")
axes[1].hist(rs_col, bins=range(1,50))
axes[1].set_xlabel("Read support (RS)")
plt.tight_layout()
plt.savefig("p_thresholds.png")