"""
TE Family Presence/Absence Analysis from VCF
=============================================
Analyses distribution of specific TE families across Z. tritici strains.
VCF uses 0/1 presence/absence genotypes per strain.

Families of interest:
  - ZymTri_2023_family_1222#LINE
  - ZymTri_2023_family_148#LINE
  - ZymTri_2023_family_164#LINE
  - ZymTri_2023_family_1058#LINE/I-Jockey
"""

# %%
# ── Imports ──────────────────────────────────────────────────────────────────
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import seaborn as sns
from pathlib import Path
from collections import defaultdict

# %%
# ── Configuration ─────────────────────────────────────────────────────────────
VCF_PATH = Path("SupplementaryData3.vcf")   # <-- adjust if needed

FAMILIES_OF_INTEREST = [
    "ZymTri_2023_family_1222#LINE",
    "ZymTri_2023_family_148#LINE",
    "ZymTri_2023_family_164#LINE",
    "ZymTri_2023_family_1058#LINE/I-Jockey",
]

# Short display labels (strips ZymTri_2023_ prefix for readability in plots)
FAMILY_LABELS = {f: f.replace("ZymTri_2023_", "") for f in FAMILIES_OF_INTEREST}

# %%
# ── VCF Parser ────────────────────────────────────────────────────────────────
def parse_vcf(vcf_path: Path):
    """
    Parse a presence/absence VCF file.

    Returns
    -------
    samples : list[str]
        Ordered list of sample/strain names.
    records : list[dict]
        One dict per variant with keys:
          'chrom', 'pos', 'id', 'ref', 'alt', 'info',
          'family'  (TE family extracted from ID or INFO),
          'genotypes' (dict sample -> 0 or 1)
    """
    samples = []
    records = []

    with open(vcf_path, "r") as fh:
        for line in fh:
            line = line.rstrip("\n")

            # ── Header: collect sample names ──────────────────────────────
            if line.startswith("#CHROM"):
                cols = line.lstrip("#").split("\t")
                # Standard VCF: CHROM POS ID REF ALT QUAL FILTER INFO FORMAT sample1 sample2 ...
                if "FORMAT" in cols:
                    fmt_idx = cols.index("FORMAT")
                    samples = cols[fmt_idx + 1:]
                else:
                    # No FORMAT column — genotypes start after INFO (col 8)
                    samples = cols[8:]
                continue

            if line.startswith("#"):
                continue  # other meta-lines

            # ── Data line ────────────────────────────────────────────────
            parts = line.split("\t")
            chrom, pos, var_id, ref, alt = parts[0], int(parts[1]), parts[2], parts[3], parts[4]
            info = parts[7] if len(parts) > 7 else ""

            # Determine family: prefer ID field, fall back to INFO SVTYPE/TE= tag
            family = _extract_family(var_id, info)

            # Parse genotypes — handle GT-only (0/1) or GT:... formats
            has_format = len(parts) > 8 and ":" in parts[8] if len(parts) > 8 else False
            gt_start = 9 if has_format else 8

            gts = {}
            for i, sample in enumerate(samples):
                raw = parts[gt_start + i] if gt_start + i < len(parts) else "."
                gts[sample] = _parse_gt(raw)

            records.append({
                "chrom": chrom,
                "pos": pos,
                "id": var_id,
                "ref": ref,
                "alt": alt,
                "info": info,
                "family": family,
                "genotypes": gts,
            })

    return samples, records


def _extract_family(var_id: str, info: str) -> str:
    """
    Extract TE family name from the variant ID or INFO field.
    Tries common patterns:
      - ID contains '#' (TEsorter-style: name#superfamily)
      - INFO has FAMILY=... or TE=... or SVTYPE=...
    """
    # ID field often encodes family directly
    if var_id and var_id != ".":
        # Strip copy-number suffix, e.g. "ZymTri_2023_family_1222#LINE_copy001"
        # Keep everything up to the last underscore+digit block if present
        # but return the family portion (before copy suffix)
        m = re.match(r"(ZymTri_\d+_family_\d+#[^_\s]+)", var_id)
        if m:
            return m.group(1)
        return var_id  # fallback: use full ID

    # Try INFO tags
    for tag in ("FAMILY", "TE", "SVTYPE"):
        m = re.search(rf"{tag}=([^;]+)", info)
        if m:
            return m.group(1)

    return "unknown"


def _parse_gt(raw: str) -> int:
    """Return 1 if any alt allele present, 0 if homRef, -1 if missing."""
    gt = raw.split(":")[0]       # take GT field only
    if "." in gt:
        return -1                # missing
    alleles = re.split(r"[/|]", gt)
    return 1 if any(a != "0" for a in alleles) else 0


# %%
# ── Load Data ─────────────────────────────────────────────────────────────────
print(f"Parsing {VCF_PATH} …")
samples, records = parse_vcf(VCF_PATH)
print(f"  Strains : {len(samples)}")
print(f"  Variants: {len(records)}")

# %%
# ── Build Tidy DataFrame ──────────────────────────────────────────────────────
rows = []
for rec in records:
    for sample, gt in rec["genotypes"].items():
        rows.append({
            "variant_id": rec["id"],
            "family":     rec["family"],
            "chrom":      rec["chrom"],
            "pos":        rec["pos"],
            "strain":     sample,
            "presence":   gt,   # 1=present, 0=absent, -1=missing
        })

df_all = pd.DataFrame(rows)
print(f"\nFull tidy dataframe shape: {df_all.shape}")
print(df_all.head())

# %%
# ── Filter to Families of Interest ───────────────────────────────────────────
df = df_all[df_all["family"].isin(FAMILIES_OF_INTEREST)].copy()
df["family_label"] = df["family"].map(FAMILY_LABELS)

print(f"\nRows for families of interest: {len(df)}")
print(df["family"].value_counts())

# %%
# ── Summary Table: copies per family ─────────────────────────────────────────
# Each unique variant_id = one TE copy
copies_per_family = (
    df.groupby("family")["variant_id"]
    .nunique()
    .rename("n_copies")
    .reset_index()
)
copies_per_family["family_label"] = copies_per_family["family"].map(FAMILY_LABELS)
print("\n── Copies (TE insertions) per family ──")
print(copies_per_family[["family_label", "n_copies"]].to_string(index=False))

# %%
# ── Summary Table: presence/absence per family × strain ──────────────────────
# For each family and strain: how many copies are present?
presence_counts = (
    df[df["presence"] == 1]
    .groupby(["family", "strain"])["variant_id"]
    .nunique()
    .reset_index(name="copies_present")
)

# Fill strains that have zero copies present
full_index = pd.MultiIndex.from_product(
    [FAMILIES_OF_INTEREST, samples], names=["family", "strain"]
)
presence_counts = (
    presence_counts
    .set_index(["family", "strain"])
    .reindex(full_index, fill_value=0)
    .reset_index()
)
presence_counts["family_label"] = presence_counts["family"].map(FAMILY_LABELS)

print("\n── Copy counts per strain (first 10 rows) ──")
print(presence_counts.head(10).to_string(index=False))

# %%
# ── Binary presence matrix (copy-level) ──────────────────────────────────────
# Rows = TE copies, Columns = strains, Values = 0/1
df_pres = df[df["presence"].isin([0, 1])].copy()
pres_matrix = df_pres.pivot_table(
    index=["family", "variant_id"],
    columns="strain",
    values="presence",
    fill_value=0,
).astype(int)

print(f"\nPresence matrix shape: {pres_matrix.shape}")
print(pres_matrix.iloc[:5, :5])

# %%
# ── Per-family statistics ─────────────────────────────────────────────────────
stats_rows = []
for fam in FAMILIES_OF_INTEREST:
    sub = df[(df["family"] == fam) & (df["presence"].isin([0, 1]))]
    copies = sub["variant_id"].unique()
    n_copies = len(copies)

    # Per-copy: fraction of strains where present
    copy_prev = sub.groupby("variant_id")["presence"].mean()  # prevalence per copy

    # Per-strain: how many copies present
    strain_counts = (
        sub[sub["presence"] == 1]
        .groupby("strain")["variant_id"]
        .nunique()
        .reindex(samples, fill_value=0)
    )

    stats_rows.append({
        "family":                   fam,
        "family_label":             FAMILY_LABELS[fam],
        "n_copies":                 n_copies,
        "n_strains_total":          len(samples),
        "n_strains_any_present":    (strain_counts > 0).sum(),
        "pct_strains_any_present":  round(100 * (strain_counts > 0).mean(), 1),
        "mean_copies_per_strain":   round(strain_counts.mean(), 2),
        "median_copies_per_strain": round(strain_counts.median(), 2),
        "max_copies_per_strain":    int(strain_counts.max()),
        "mean_copy_prevalence":     round(copy_prev.mean(), 3),
    })

stats_df = pd.DataFrame(stats_rows)
print("\n── Per-family statistics ──")
print(stats_df.to_string(index=False))
stats_df.to_csv("te_family_stats.csv", index=False)
print("Saved: te_family_stats.csv")

# %%
# ── Plot 1: Copy counts per family ───────────────────────────────────────────
fig, ax = plt.subplots(figsize=(7, 4))
colors = sns.color_palette("Set2", len(FAMILIES_OF_INTEREST))
bars = ax.bar(
    stats_df["family_label"],
    stats_df["n_copies"],
    color=colors,
    edgecolor="black",
    linewidth=0.7,
)
ax.bar_label(bars, padding=3, fontsize=10)
ax.set_title("Number of TE copies per family", fontsize=13, fontweight="bold")
ax.set_xlabel("Family", fontsize=11)
ax.set_ylabel("Number of copies (insertions)", fontsize=11)
ax.tick_params(axis="x", rotation=20)
plt.tight_layout()
plt.savefig("plot1_copies_per_family.pdf", dpi=150)
plt.savefig("plot1_copies_per_family.png", dpi=150)
plt.show()
print("Saved: plot1_copies_per_family.pdf / .png")

# %%
# ── Plot 2: Distribution of copy counts across strains (per family) ───────────
fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharey=False)
axes = axes.flatten()

for ax, fam, color in zip(axes, FAMILIES_OF_INTEREST, colors):
    sub = df[(df["family"] == fam) & (df["presence"].isin([0, 1]))]
    strain_counts = (
        sub[sub["presence"] == 1]
        .groupby("strain")["variant_id"]
        .nunique()
        .reindex(samples, fill_value=0)
    )
    ax.hist(strain_counts, bins=range(0, strain_counts.max() + 2),
            color=color, edgecolor="black", linewidth=0.6, align="left")
    ax.set_title(FAMILY_LABELS[fam], fontsize=10, fontweight="bold")
    ax.set_xlabel("Copies present per strain", fontsize=9)
    ax.set_ylabel("Number of strains", fontsize=9)
    ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))
    # Annotate median
    med = strain_counts.median()
    ax.axvline(med, color="red", linestyle="--", linewidth=1.2, label=f"Median={med:.1f}")
    ax.legend(fontsize=8)

plt.suptitle("Distribution of TE copy counts across strains", fontsize=13, fontweight="bold", y=1.01)
plt.tight_layout()
plt.savefig("plot2_copy_count_distributions.pdf", dpi=150, bbox_inches="tight")
plt.savefig("plot2_copy_count_distributions.png", dpi=150, bbox_inches="tight")
plt.show()
print("Saved: plot2_copy_count_distributions.pdf / .png")

# %%
# ── Plot 3: Heatmap – copy-level presence/absence across strains ──────────────
# One panel per family (too many copies to combine)
for fam in FAMILIES_OF_INTEREST:
    sub_matrix = pres_matrix.loc[fam] if fam in pres_matrix.index.get_level_values("family") else None
    if sub_matrix is None or sub_matrix.empty:
        print(f"No data for {fam}, skipping heatmap.")
        continue

    n_copies, n_strains = sub_matrix.shape
    fig_h = max(4, n_copies * 0.25)
    fig_w = max(8, n_strains * 0.15)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    sns.heatmap(
        sub_matrix,
        ax=ax,
        cmap=["#f0f0f0", "#2d6a4f"],
        cbar_kws={"label": "Presence (1) / Absence (0)", "shrink": 0.5},
        linewidths=0 if n_copies > 50 else 0.3,
        yticklabels=n_copies <= 80,
        xticklabels=n_strains <= 80,
        vmin=0, vmax=1,
    )
    short = FAMILY_LABELS[fam]
    ax.set_title(f"Presence/absence: {short}", fontsize=11, fontweight="bold")
    ax.set_xlabel("Strain", fontsize=9)
    ax.set_ylabel("TE copy", fontsize=9)
    plt.tight_layout()

    safe_name = re.sub(r"[^\w]", "_", short)
    plt.savefig(f"plot3_heatmap_{safe_name}.pdf", dpi=150, bbox_inches="tight")
    plt.savefig(f"plot3_heatmap_{safe_name}.png", dpi=150, bbox_inches="tight")
    plt.show()
    print(f"Saved: plot3_heatmap_{safe_name}.pdf / .png")

# %%
# ── Plot 4: Prevalence per copy (fraction of strains where present) ───────────
fig, ax = plt.subplots(figsize=(7, 4))
for fam, color in zip(FAMILIES_OF_INTEREST, colors):
    sub = df[(df["family"] == fam) & (df["presence"].isin([0, 1]))]
    copy_prev = sub.groupby("variant_id")["presence"].mean().sort_values()
    ax.plot(
        range(len(copy_prev)),
        copy_prev.values,
        marker="o", markersize=3, linewidth=1,
        label=FAMILY_LABELS[fam], color=color,
    )

ax.set_title("Copy prevalence (fraction of strains where present)", fontsize=12, fontweight="bold")
ax.set_xlabel("TE copy rank", fontsize=10)
ax.set_ylabel("Fraction of strains with copy", fontsize=10)
ax.set_ylim(0, 1.05)
ax.legend(fontsize=8, loc="upper left")
plt.tight_layout()
plt.savefig("plot4_copy_prevalence.pdf", dpi=150)
plt.savefig("plot4_copy_prevalence.png", dpi=150)
plt.show()
print("Saved: plot4_copy_prevalence.pdf / .png")

# %%
# ── Plot 5: Strain-level presence across families (stacked bar) ───────────────
# For each strain, how many copies of each family are present?
strain_family_counts = (
    df[df["presence"] == 1]
    .groupby(["strain", "family"])["variant_id"]
    .nunique()
    .unstack(fill_value=0)
    .reindex(columns=FAMILIES_OF_INTEREST, fill_value=0)
    .reindex(samples, fill_value=0)
)
strain_family_counts.columns = [FAMILY_LABELS[c] for c in strain_family_counts.columns]

# Sort strains by total TE burden
strain_family_counts = strain_family_counts.loc[
    strain_family_counts.sum(axis=1).sort_values(ascending=False).index
]

fig, ax = plt.subplots(figsize=(max(10, len(samples) * 0.18), 5))
strain_family_counts.plot(
    kind="bar", stacked=True, ax=ax,
    color=sns.color_palette("Set2", len(FAMILIES_OF_INTEREST)),
    edgecolor="none", width=0.85,
)
ax.set_title("TE copy burden per strain (families of interest)", fontsize=12, fontweight="bold")
ax.set_xlabel("Strain", fontsize=10)
ax.set_ylabel("Number of TE copies present", fontsize=10)
ax.legend(title="Family", bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8)
ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=6)
plt.tight_layout()
plt.savefig("plot5_per_strain_stacked.pdf", dpi=150, bbox_inches="tight")
plt.savefig("plot5_per_strain_stacked.png", dpi=150, bbox_inches="tight")
plt.show()
print("Saved: plot5_per_strain_stacked.pdf / .png")

# %%
# ── Export tidy presence table ────────────────────────────────────────────────
out = (
    df[df["presence"].isin([0, 1])]
    [["family", "family_label", "variant_id", "chrom", "pos", "strain", "presence"]]
    .sort_values(["family", "variant_id", "strain"])
)
out.to_csv("te_presence_tidy.csv", index=False)
print("Saved: te_presence_tidy.csv")

pres_matrix.to_csv("te_presence_matrix.csv")
print("Saved: te_presence_matrix.csv  (copies × strains binary matrix)")

# %%
print("\n✅ All done. Output files:")
print("  te_family_stats.csv")
print("  te_presence_tidy.csv")
print("  te_presence_matrix.csv")
print("  plot1_copies_per_family.{pdf,png}")
print("  plot2_copy_count_distributions.{pdf,png}")
print("  plot3_heatmap_<family>.{pdf,png}  (one per family)")
print("  plot4_copy_prevalence.{pdf,png}")
print("  plot5_per_strain_stacked.{pdf,png}")