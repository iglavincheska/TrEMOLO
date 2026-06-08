#%%
import pandas as pd
import numpy as np
from pathlib import Path
import sys

# Define your sample directories here
SAMPLES = {
    #"IPO323_WT_18C_R1_sc": "/Users/iglavinchesk/packages/TrEMOLO/IPO323/results/IPO323_WT_18C_R1_sc_vs_parent_only",
    #"IPO323_WT_18C_R2_sc": "/Users/iglavinchesk/packages/TrEMOLO/IPO323/results/IPO323_WT_18C_R2_sc_vs_parent_only",
    #"IPO323_WT_18C_R3_sc": "/Users/iglavinchesk/packages/TrEMOLO/IPO323/results/IPO323_WT_18C_R3_sc_vs_parent_only",
    "IPO323_WT_18C_R1_p": "/Users/iglavinchesk/packages/TrEMOLO/IPO323/results/IPO323_WT_18C_R1_p_vs_parent_only",
    "IPO323_WT_18C_R2_p": "/Users/iglavinchesk/packages/TrEMOLO/IPO323/results/IPO323_WT_18C_R2_p_vs_parent_only",
    "IPO323_WT_18C_R3_p": "/Users/iglavinchesk/packages/TrEMOLO/IPO323/results/IPO323_WT_18C_R3_p_vs_parent_only",
    
    "IPO323_dim2_18C_R1_p": "/Users/iglavinchesk/packages/TrEMOLO/IPO323/results/IPO323_dim2_18C_R1_p_vs_parent_only",
    "IPO323_dim2_18C_R3_p": "/Users/iglavinchesk/packages/TrEMOLO/IPO323/results/IPO323_dim2_18C_R3_p_vs_parent_only",
    
    "IPO323_WT_28C_R1_p": "/Users/iglavinchesk/packages/TrEMOLO/IPO323/results/IPO323_WT_28C_R1_p_vs_parent_only",
    "IPO323_WT_28C_R2_p": "/Users/iglavinchesk/packages/TrEMOLO/IPO323/results/IPO323_WT_28C_R2_p_vs_parent_only",
    "IPO323_WT_28C_R3_p": "/Users/iglavinchesk/packages/TrEMOLO/IPO323/results/IPO323_WT_28C_R3_p_vs_parent_only",
    "IPO323_dim2_28C_R1_p": "/Users/iglavinchesk/packages/TrEMOLO/IPO323/results/IPO323_dim2_28C_R1_p_vs_parent_only",
    "IPO323_dim2_28C_R2_p": "/Users/iglavinchesk/packages/TrEMOLO/IPO323/results/IPO323_dim2_28C_R2_p_vs_parent_only",
    "IPO323_dim2_28C_R3_p": "/Users/iglavinchesk/packages/TrEMOLO/IPO323/results/IPO323_dim2_28C_R3_p_vs_parent_only",
    
    #"3D7_WT_18C_R1_p": "/Users/iglavinchesk/packages/TrEMOLO/3D7/results/3D7_WT_18C_R1_p_vs_parent_only",
    "3D7_WT_18C_R2_p": "/Users/iglavinchesk/packages/TrEMOLO/3D7/results/3D7_WT_18C_R2_p_vs_parent_only",
    "3D7_WT_18C_R3_p": "/Users/iglavinchesk/packages/TrEMOLO/3D7/results/3D7_WT_18C_R3_p_vs_parent_only",
    
    #"3D7_dim2_18C_R1_p": "/Users/iglavinchesk/packages/TrEMOLO/3D7/results/3D7_dim2_18C_R1_p_vs_parent_only",
    "3D7_dim2_18C_R2_p": "/Users/iglavinchesk/packages/TrEMOLO/3D7/results/3D7_dim2_18C_R2_p_vs_parent_only",
    "3D7_dim2_18C_R3_p": "/Users/iglavinchesk/packages/TrEMOLO/3D7/results/3D7_dim2_18C_R3_p_vs_parent_only",
    
    "IR0148B_WT_18C_R1_p": "/Users/iglavinchesk/packages/TrEMOLO/IR0148B/results/IR0148B_WT_18C_R1_p_vs_parent_only",
    "IR0148B_WT_18C_R2_p": "/Users/iglavinchesk/packages/TrEMOLO/IR0148B/results/IR0148B_WT_18C_R2_p_vs_parent_only",
    "IR0148B_WT_18C_R3_p": "/Users/iglavinchesk/packages/TrEMOLO/IR0148B/results/IR0148B_WT_18C_R3_p_vs_parent_only",
    
    "IR0126B_WT_18C_R1_p": "/Users/iglavinchesk/packages/TrEMOLO/IR0126B/results/IR0126B_WT_18C_R1_p_vs_parent_only",
    "IR0126B_WT_18C_R2_p": "/Users/iglavinchesk/packages/TrEMOLO/IR0126B/results/IR0126B_WT_18C_R2_p_vs_parent_only",
    "IR0126B_WT_18C_R3_p": "/Users/iglavinchesk/packages/TrEMOLO/IR0126B/results/IR0126B_WT_18C_R3_p_vs_parent_only",
    "IR0126B_WT_18C_R1_sc": "/Users/iglavinchesk/packages/TrEMOLO/IR0126B/results/IR0126B_WT_18C_R1_sc_vs_ref_insider",
    "IR0126B_WT_18C_R2_sc": "/Users/iglavinchesk/packages/TrEMOLO/IR0126B/results/IR0126B_WT_18C_R2_sc_vs_ref_insider",
    
    "IR0126B_dim2_18C_R1_p": "/Users/iglavinchesk/packages/TrEMOLO/IR0126B/results/IR0126B_dim2_18C_R1_p_vs_parent_only",
    "IR0126B_dim2_18C_R2_p": "/Users/iglavinchesk/packages/TrEMOLO/IR0126B/results/IR0126B_dim2_18C_R2_p_vs_parent_only",
    "IR0126B_dim2_18C_R3_p": "/Users/iglavinchesk/packages/TrEMOLO/IR0126B/results/IR0126B_dim2_18C_R3_p_vs_parent_only",
}


def load_te_bed(filepath):
    """Load TE_INFOS.bed and return as DataFrame."""
    try:
        df = pd.read_csv(filepath, sep="\t", header=0)
        return df
    except FileNotFoundError:
        print(f"  ⚠ File not found: {filepath}")
        return None


def filter_outsider(df):
    """Keep only OUTSIDER TE insertions."""
    return df[~df["ID_TrEMOLO"].str.contains("INSIDER", na=False)].copy()


def summarize_sample(sample_name, df):
    """Generate summary statistics for a single sample."""
    if df is None:
        return None

    # Include even if empty (0 TEs)
    outsider = filter_outsider(df)
    insider = df[df["ID_TrEMOLO"].str.contains("INSIDER", na=False)] if not df.empty else pd.DataFrame()

    summary = {
        "Sample": sample_name,
        "Total_TEs": len(df),
        "Outsider_TEs": len(outsider),
        "Insider_TEs": len(insider),
    }

    # Add OUTSIDER statistics only if there are OUTSIDER TEs
    if not outsider.empty:
        freq_col = pd.to_numeric(outsider["FREQ_WITH_CLIPPED"], errors="coerce")
        if "RS" in outsider.columns:
            rs_col = pd.to_numeric(outsider["RS"], errors="coerce")
        else:
            rs_col = pd.Series([np.nan] * len(outsider))

        summary.update({
            "Outsider_FREQ≥5%": (freq_col >= 5).sum(),
            "Outsider_FREQ≥10%": (freq_col >= 10).sum(),
            "Outsider_FREQ≥20%": (freq_col >= 20).sum(),
            "Outsider_RS≥2": (rs_col >= 2).sum(),
            "Outsider_RS≥3": (rs_col >= 3).sum(),
            "Outsider_RS≥5": (rs_col >= 5).sum(),
            "FREQ_mean": freq_col.mean(),
            "FREQ_median": freq_col.median(),
            "FREQ_std": freq_col.std(),
            "FREQ_max": freq_col.max(),
            "RS_mean": rs_col.mean(),
            "RS_median": rs_col.median(),
            "RS_max": rs_col.max(),
        })
    else:
        # Fill with 0 for samples with no OUTSIDER TEs
        summary.update({
            "Outsider_FREQ≥5%": 0,
            "Outsider_FREQ≥10%": 0,
            "Outsider_FREQ≥20%": 0,
            "Outsider_RS≥2": 0,
            "Outsider_RS≥3": 0,
            "Outsider_RS≥5": 0,
            "FREQ_mean": np.nan,
            "FREQ_median": np.nan,
            "FREQ_std": np.nan,
            "FREQ_max": np.nan,
            "RS_mean": np.nan,
            "RS_median": np.nan,
            "RS_max": np.nan,
        })

    return summary

#%%

print("=" * 80)
print("TrEMOLO TE_INFOS.bed Analysis")
print("=" * 80)

results = []
all_data = {}

    # Load and analyze each sample
for sample_name, sample_dir in SAMPLES.items():
    te_bed_path = Path(sample_dir) / "TE_INFOS.bed"
    print(f"\n📊 Loading {sample_name}...")

    df = load_te_bed(str(te_bed_path))
    if df is not None:
        all_data[sample_name] = filter_outsider(df)
        summary = summarize_sample(sample_name, df)
        if summary:
            results.append(summary)
            insider = df[df["ID_TrEMOLO"].str.contains("INSIDER", na=False)] if not df.empty else pd.DataFrame()
            print(f"   ✓ Loaded {len(df)} TEs (outsider: {len(all_data[sample_name])}, insider: {len(insider)})")

if not results:
    print("\n❌ No valid samples loaded.")
    sys.exit(1)

    # Create summary table
summary_df = pd.DataFrame(results)
print("\n" + "=" * 80)
print("SUMMARY TABLE")
print("=" * 80)
print(summary_df.to_string(index=False))

    # Save summary table
summary_df.to_csv("/Users/iglavinchesk/packages/TrEMOLO/IPO323/results/tremolo_summary.tsv", sep="\t", index=False)
print(f"\n✓ Summary saved to /Users/iglavinchesk/packages/TrEMOLO/IPO323/results/tremolo_summary.tsv")


#%%

    # Comparison: Population (p) vs Single Colony (sc)
print("\n" + "=" * 80)
print("POPULATION vs SINGLE COLONY COMPARISON")
print("=" * 80)

comparisons = [
    ("R1", "IPO323_WT_18C_R1_p", "IPO323_WT_18C_R1_sc"),
    ("R2", "IPO323_WT_18C_R2_p", "IPO323_WT_18C_R2_sc"),
    ("R3", "IPO323_WT_18C_R3_p", "IPO323_WT_18C_R3_sc"),
    ]

for rep, pop_name, sc_name in comparisons:
    if pop_name in all_data and sc_name in all_data:
        pop_df = all_data[pop_name]
        sc_df = all_data[sc_name]

        print(f"\n{rep}:")
        print(f"  Population:")
        print(f"    - Total outsider TEs: {len(pop_df)}")
        print(f"    - FREQ ≥ 5%: {(pd.to_numeric(pop_df['FREQ_WITH_CLIPPED'], errors='coerce') >= 5).sum()}")
        if "RS" in pop_df.columns:
            print(f"    - RS ≥ 2: {(pd.to_numeric(pop_df['RS'], errors='coerce') >= 2).sum()}")
            print(f"    - RS ≥ 3: {(pd.to_numeric(pop_df['RS'], errors='coerce') >= 3).sum()}")

        print(f"  Single Colony:")
        print(f"    - Total outsider TEs: {len(sc_df)}")
        print(f"    - FREQ ≥ 5%: {(pd.to_numeric(sc_df['FREQ_WITH_CLIPPED'], errors='coerce') >= 5).sum()}")
        if "RS" in sc_df.columns:
            print(f"    - RS ≥ 2: {(pd.to_numeric(sc_df['RS'], errors='coerce') >= 2).sum()}")
            print(f"    - RS ≥ 3: {(pd.to_numeric(sc_df['RS'], errors='coerce') >= 3).sum()}")

    # Distribution of RS by sample type
has_rs = any("RS" in all_data[k].columns for k in all_data)
if has_rs:
    print("\n" + "=" * 80)
    print("RS DISTRIBUTION BY SAMPLE TYPE")
    print("=" * 80)

    pop_rs_list = []
    sc_rs_list = []
    for k in all_data:
        if "RS" in all_data[k].columns:
            rs = pd.to_numeric(all_data[k]["RS"], errors="coerce")
            if "_p" in k:
                pop_rs_list.append(rs)
            elif "_sc" in k:
                sc_rs_list.append(rs)

    if pop_rs_list:
        pop_rs = pd.concat(pop_rs_list, ignore_index=True)
        print("\nPopulation samples (p):")
        print(pop_rs.describe())

    if sc_rs_list:
        sc_rs = pd.concat(sc_rs_list, ignore_index=True)
        print("\nSingle Colony samples (sc):")
        print(sc_rs.describe())
else:
    print("\n⚠ RS column not available in loaded files (old format)")

#%%

    # Frequency thresholds summary
print("\n" + "=" * 80)
print("SUGGESTED FILTERING BY THRESHOLD")
print("=" * 80)

thresholds = {
    "MIN_FREQ": [0, 5, 10, 20],
    "MIN_RS": [1, 2, 3, 5],
}

for sample_name, df in list(all_data.items()):  # Show all samples
    if df.empty:
        continue  # Skip samples with no TEs

    print(f"\n{sample_name}:")
    outsider = filter_outsider(df)

    if outsider.empty:
        print("  No OUTSIDER TEs")
        continue

    freq = pd.to_numeric(outsider["FREQ_WITH_CLIPPED"], errors="coerce")

    print("  MIN_FREQ thresholds (% frequency):")
    for threshold in thresholds["MIN_FREQ"]:
        count = (freq >= threshold).sum()
        pct = 100 * count / len(outsider)
        print(f"    ≥ {threshold:2d}%: {count:3d} TEs ({pct:5.1f}%)")

    if "RS" in df.columns:
        rs = pd.to_numeric(outsider["RS"], errors="coerce")
        print("  MIN_RS thresholds (read support):")
        for threshold in thresholds["MIN_RS"]:
            count = (rs >= threshold).sum()
            pct = 100 * count / len(outsider)
            print(f"    ≥ {threshold}: {count:3d} TEs ({pct:5.1f}%)")
    else:
        print("  MIN_RS: (RS column not available in this file)")

#%%

# High-confidence TEs (RS >= 2 AND FREQ > 5%) - Separated by type
print("\n" + "=" * 80)
print("HIGH-CONFIDENCE TEs (RS ≥ 2 AND FREQ > 5%) - By Type")
print("=" * 80)

import matplotlib.pyplot as plt

hc_data = []
hc_by_strain = {}  # {strain: {ins: count, del: count, other: count}}

# First pass: collect high-confidence TEs and build hc_by_strain
hc_samples = {}  # Track which samples have high-confidence data
for sample_name, df in all_data.items():
    if df.empty:
        continue

    outsider = filter_outsider(df)

    # Extract strain name
    strain = sample_name.split("_")[0]

    if strain not in hc_by_strain:
        hc_by_strain[strain] = {"INS": 0, "DEL": 0, "OTHER": 0}

    if not outsider.empty and "RS" in outsider.columns:
        freq_col = pd.to_numeric(outsider["FREQ_WITH_CLIPPED"], errors="coerce")
        rs_col = pd.to_numeric(outsider["RS"], errors="coerce")
        hc = outsider[(rs_col >= 2) & (freq_col > 5)]

        # Count by TYPE
        ins_count = len(hc[hc["TYPE"].str.contains("INS", na=False)]) if "TYPE" in hc.columns else 0
        del_count = len(hc[hc["TYPE"].str.contains("DEL", na=False)]) if "TYPE" in hc.columns else 0
        other_count = len(hc) - ins_count - del_count

        hc_by_strain[strain]["INS"] += ins_count
        hc_by_strain[strain]["DEL"] += del_count
        hc_by_strain[strain]["OTHER"] += other_count

        hc_data.append({
            "Sample": sample_name,
            "Strain": strain,
            "INS": ins_count,
            "DEL": del_count,
            "OTHER": other_count,
            "Total": len(hc)
        })
        hc_samples[sample_name] = True

        print(f"\n{sample_name}:")
        print(f"  INS: {ins_count}, DEL: {del_count}, OTHER: {other_count}, TOTAL: {len(hc)}")

# Second pass: add zero-value entries for samples without high-confidence TEs
for sample_name, df in all_data.items():
    if sample_name not in hc_samples:
        strain = sample_name.split("_")[0]
        hc_data.append({
            "Sample": sample_name,
            "Strain": strain,
            "INS": 0,
            "DEL": 0,
            "OTHER": 0,
            "Total": 0
        })

# Create visualization - organized by STRAIN, with replicates as adjacent bars
if hc_data:
    hc_df = pd.DataFrame(hc_data)

    # Parse sample info: Strain_Type_Temp_Rep_Condition
    hc_df["Strain"] = hc_df["Sample"].apply(lambda x: x.split("_")[0])
    hc_df["Type"] = hc_df["Sample"].apply(lambda x: "dim2" if "dim2" in x else "WT")
    hc_df["Replicate"] = hc_df["Sample"].apply(lambda x: x.split("_")[-2] if "_" in x else "")
    hc_df["Temperature"] = hc_df["Sample"].apply(lambda x: [s for s in x.split("_") if "C" in s][0] if any("C" in s for s in x.split("_")) else "")

    # Get unique values and sort them
    replicates = sorted(hc_df["Replicate"].unique())
    strains = sorted(hc_df["Strain"].unique())
    types = ["WT", "dim2"]

    # Create a single plot with strains on X-axis and replicates as adjacent bars
    fig, ax = plt.subplots(figsize=(16, 6))

    # Colors and hatches
    colors_type = {"WT": "steelblue", "dim2": "darkorange"}
    hatch_del = "///"
    hatch_other = "xxx"

    # Calculate bar positions
    bar_width = 0.12  # Width of each individual bar
    rep_spacing = 0.05  # Space between WT and dim2 within replicate
    group_spacing = 0.5  # Space between replicates within strain
    strain_spacing = 0.8  # Space between strain groups

    # Position for each (strain, replicate, type) combination
    x_pos = 0
    strain_x_positions = {}  # Track center position of each strain group

    for strain_idx, strain in enumerate(strains):
        strain_data = hc_df[hc_df["Strain"] == strain]
        if strain_data.empty:
            continue

        strain_start = x_pos
        strain_positions = []

        # For each replicate
        for rep in replicates:
            rep_data = strain_data[strain_data["Replicate"] == rep]
            if rep_data.empty:
                continue

            rep_start = x_pos

            # For each type (WT, dim2) - always show both side by side
            for type_val in types:
                type_data = rep_data[rep_data["Type"] == type_val]

                # Check if data exists for this combination
                if not type_data.empty:
                    row = type_data.iloc[0]
                    is_zero = row["Total"] == 0
                else:
                    # Create a zero-value row for missing data
                    row = {"INS": 0, "DEL": 0, "OTHER": 0, "Total": 0}
                    is_zero = True

                # Use lighter/different color for zero values
                if is_zero:
                    bar_color = colors_type[type_val]
                    bar_alpha = 0.3
                    # Show a small indicator bar even for zeros
                    min_height = 0.3
                else:
                    bar_color = colors_type[type_val]
                    bar_alpha = 1.0
                    min_height = 0

                # Stacked bars: INS, DEL, OTHER
                total_height = max(row["INS"] + row["DEL"] + row["OTHER"], min_height)

                ax.bar(x_pos, row["INS"], bar_width, color=bar_color, alpha=bar_alpha,
                       edgecolor="black", linewidth=0.5)
                ax.bar(x_pos, row["DEL"], bar_width, bottom=row["INS"],
                       color=bar_color, hatch=hatch_del, alpha=bar_alpha,
                       edgecolor="black", linewidth=0.5)
                ax.bar(x_pos, row["OTHER"], bar_width, bottom=row["INS"] + row["DEL"],
                       color=bar_color, hatch=hatch_other, alpha=bar_alpha,
                       edgecolor="black", linewidth=0.5)

                # Add total on top (show "0" for zero values)
                total = row["Total"]
                label_y = total_height + 0.5 if total > 0 else min_height + 0.3
                ax.text(x_pos, label_y, str(int(total)), ha="center", va="bottom",
                       fontweight="bold", fontsize=8)

                strain_positions.append(x_pos)
                x_pos += bar_width

                # Small gap between WT and dim2
                if type_val == "WT":
                    x_pos += rep_spacing

            # Gap between replicates within strain
            x_pos += group_spacing

        # Record center position for strain label
        if strain_positions:
            strain_x_positions[strain] = (min(strain_positions) + max(strain_positions)) / 2

        # Gap between strains
        x_pos += strain_spacing

    # Set X-axis labels at strain level
    ax.set_xticks([pos for pos in strain_x_positions.values()])
    ax.set_xticklabels([strain for strain in strains], fontsize=11, fontweight="bold")
    ax.set_ylabel("Number of TEs", fontsize=12, fontweight="bold")
    ax.set_title("High-Confidence TEs by Strain and Replicate\n(RS ≥ 2, FREQ > 5%)",
                fontsize=13, fontweight="bold")
    ax.grid(axis="y", alpha=0.3)
    ax.set_ylim(0, max(hc_df["Total"]) * 1.15)

    # Custom legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="steelblue", edgecolor="black", label="WT"),
        Patch(facecolor="darkorange", edgecolor="black", label="dim2"),
        Patch(facecolor="gray", edgecolor="black", label="INS (solid)"),
        Patch(facecolor="gray", hatch="///", edgecolor="black", label="DEL (hatched)"),
        Patch(facecolor="gray", hatch="xxx", edgecolor="black", alpha=0.6, label="OTHER (cross-hatched)"),
    ]
    ax.legend(handles=legend_elements, loc="upper right", fontsize=10, ncol=2)

    plt.tight_layout()
    plt.savefig("/Users/iglavinchesk/packages/TrEMOLO/downstream_analysis/hc_tes_by_strain.png", dpi=300, bbox_inches="tight")
    print("\n✓ Plot saved to /Users/iglavinchesk/packages/TrEMOLO/downstream_analysis/hc_tes_by_strain.png")

    print("\n" + "=" * 80)
    print("DETAILED SUMMARY")
    print("=" * 80)
    print(hc_df[["Sample", "Strain", "Type", "Replicate", "INS", "DEL", "OTHER", "Total"]].to_string(index=False))
else:
    print("\n⚠ No high-confidence TEs found")

#%%

# ANALYZE CO-OCCURRING TEs (Multiple TEs in single insertions)
print("\n" + "=" * 80)
print("CO-OCCURRING TEs ANALYSIS - Investigating Multi-TE Insertions")
print("=" * 80)

# Map sample names to their result directories
sample_dirs = {}
for sample_name, sample_dir in SAMPLES.items():
    result_dir = Path(sample_dir)
    sample_dirs[sample_name] = result_dir

# Collect multi-TE insertions
multi_te_insertions = {}  # {sample: {insertion_id: [list of TEs]}}
multi_te_details = {}     # {sample: {insertion_id: [{te: ..., pident: ..., coverage: ...}, ...]}}

for sample_name, result_dir in sample_dirs.items():
    multiple_te_file = result_dir / "MULTIPLE_TE_BY_ID.txt"

    if not multiple_te_file.exists():
        continue

    try:
        df = pd.read_csv(multiple_te_file, sep="\t", header=None, names=["ins_id", "te_id", "count"])

        if df.empty:
            continue

        # Group by insertion ID
        grouped = df.groupby("ins_id")["te_id"].apply(list)

        multi_te_insertions[sample_name] = {}
        for ins_id, te_list in grouped.items():
            if len(te_list) > 1:  # Only multi-TE insertions
                multi_te_insertions[sample_name][ins_id] = te_list
    except Exception as e:
        print(f"  ⚠ Error reading {sample_name}: {e}")

if multi_te_insertions:
    print("\n" + "=" * 80)
    print("INSERTIONS WITH MULTIPLE TE HITS")
    print("=" * 80)

    # Identify most common TE combinations
    te_combinations = {}  # {frozenset(tes): count}

    for sample_name, ins_dict in multi_te_insertions.items():
        print(f"\n{sample_name}: {len(ins_dict)} insertions with multiple TE hits")

        for ins_id, te_list in sorted(ins_dict.items()):
            print(f"  {ins_id}: {', '.join(sorted(te_list))}")

            # Track combinations
            combo = frozenset(te_list)
            te_combinations[combo] = te_combinations.get(combo, 0) + 1

    print("\n" + "=" * 80)
    print("MOST COMMON TE COMBINATIONS (appearing together)")
    print("=" * 80)

    sorted_combos = sorted(te_combinations.items(), key=lambda x: x[1], reverse=True)
    for combo, count in sorted_combos:
        te_list = sorted(list(combo))
        print(f"  {' + '.join(te_list)}: {count} insertions")

    # Now look at BLAST details for these multi-TE insertions
    print("\n" + "=" * 80)
    print("BLAST ALIGNMENT DETAILS FOR MULTI-TE INSERTIONS")
    print("=" * 80)

    for sample_name, ins_dict in multi_te_insertions.items():
        if len(ins_dict) == 0:
            continue

        # Try to find the BLAST file
        result_dir = sample_dirs[sample_name]

        # Check both TrEMOLO and TE_DETECTION BLAST files
        blast_files = [
            result_dir / "OUTSIDER/TrEMOLO_SV_TE/INS/SV_INS_CLUST.bln",
            result_dir / "OUTSIDER/TE_DETECTION/BLAST_TE.bln"
        ]

        blast_file = None
        for bf in blast_files:
            if bf.exists():
                blast_file = bf
                break

        if not blast_file:
            print(f"\n  ⚠ No BLAST file found for {sample_name}")
            continue

        print(f"\n{sample_name}:")

        # Read BLAST file
        try:
            blast_df = pd.read_csv(blast_file, sep="\t", header=None)
            blast_df.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                               "qstart", "qend", "sstart", "send", "evalue", "bitscore"]

            # For each multi-TE insertion, show the BLAST hits
            for ins_id in sorted(ins_dict.keys()):
                print(f"\n  {ins_id}:")

                # Extract the insertion name from qseqid (format: chrom:...:InsID:...:...)
                matching_hits = []
                for idx, row in blast_df.iterrows():
                    qseqid = row["qseqid"]
                    if ins_id in qseqid:
                        matching_hits.append(row)

                if matching_hits:
                    for hit in matching_hits:
                        pident = hit["pident"]
                        length = hit["length"]
                        sseqid = hit["sseqid"]
                        bitscore = hit["bitscore"]
                        qstart = int(hit["qstart"])
                        qend = int(hit["qend"])
                        coverage = abs(qend - qstart)

                        print(f"    {sseqid:20s} | Pident: {pident:6.2f}% | Coverage: {coverage:4d}bp | BitScore: {bitscore:7.1f}")
        except Exception as e:
            print(f"  ⚠ Error reading BLAST for {sample_name}: {e}")

    # Summary: Are certain TE pairs/triplets inserting together consistently?
    print("\n" + "=" * 80)
    print("HYPOTHESIS: Co-insertion patterns")
    print("=" * 80)

    for combo, count in sorted_combos[:5]:  # Top 5 combinations
        te_list = sorted(list(combo))
        if count >= 2:
            print(f"\n✓ {' + '.join(te_list)} appears together in {count} insertions")
            print(f"  → Suggests possible co-transposition or nested insertion")

            # Find which samples have this combination
            samples_with_combo = []
            for sample_name, ins_dict in multi_te_insertions.items():
                for ins_id, tes in ins_dict.items():
                    if frozenset(tes) == combo:
                        samples_with_combo.append(sample_name)

            print(f"  → Found in: {', '.join(set(samples_with_combo))}")
else:
    print("\n⚠ No insertions with multiple TE hits found")
