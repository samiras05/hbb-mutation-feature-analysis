import os
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "Data"

OUT_DIR = "integrated_figures"
os.makedirs(OUT_DIR, exist_ok=True)

plt.rcParams["figure.dpi"] = 120
plt.rcParams["axes.titlesize"] = 11
plt.rcParams["axes.labelsize"] = 9
plt.rcParams["xtick.labelsize"] = 8
plt.rcParams["ytick.labelsize"] = 8


def build_variant_from_dyn(row):
    return f"{row['ref_aa']}{int(row['position_without_met_']) + 1}{row['mut_aa']}"


def build_variant_from_foldx(row):
    return f"{row['wt_aa']}{int(row['position']) + 1}{row['mut_aa']}"


def main() -> None:
    aa_path = DATA_DIR / "AA.xlsx"
    dyn_path = DATA_DIR / "DynaMut.xlsx"
    foldx_path = DATA_DIR / "FoldX.xlsx"
    merged_path = DATA_DIR / "MergedData.csv"

    aa = pd.read_excel(aa_path)
    dyn = pd.read_excel(dyn_path)
    foldx = pd.read_excel(foldx_path)
    merged = pd.read_csv(merged_path)

    dyn["variant"] = dyn.apply(build_variant_from_dyn, axis=1)
    foldx["variant"] = foldx.apply(build_variant_from_foldx, axis=1)

    merged_sub = merged[["variant", "rawscore", "phred", "revel_revel"]].drop_duplicates()

    final = aa.copy()

    dyn_cols = ["variant", "delta_vibrational_entropy", "ddg_prediction", "delta_stability_encom", "mcsm", "sdm", "duet"]
    final = final.merge(dyn[dyn_cols], on="variant", how="left")

    foldx_energy_cols = [
        "variant",
        "SD",
        "total energy",
        "Backbone Hbond",
        "Sidechain Hbond",
        "Van der Waals",
        "Electrostatics",
        "Solvation Polar",
        "Solvation Hydrophobic",
        "Van der Waals clashes",
        "entropy sidechain",
        "entropy mainchain",
        "torsional clash",
        "backbone clash",
        "helix dipole",
        "electrostatic kon",
        "energy Ionisation",
    ]
    foldx_energy_cols = [c for c in foldx_energy_cols if c in foldx.columns]
    final = final.merge(foldx[foldx_energy_cols], on="variant", how="left")
    final = final.merge(merged_sub, on="variant", how="left")

    print("Integrated df shape:", final.shape)

    aa_variants = list(aa["variant"].astype(str))

    positions_plot = []
    energies_plot = []

    for v in aa_variants:
        m = re.match(r"^([A-Z])(\d+)([A-Z])$", v)
        if not m:
            continue
        wt, pos_aa, mut = m.group(1), int(m.group(2)), m.group(3)
        pos_fx = pos_aa - 1
        match_rows = foldx[
            (foldx["position"] == pos_fx) &
            (foldx["wt_aa"] == wt) &
            (foldx["mut_aa"] == mut)
        ]
        if len(match_rows) == 1 and "total energy" in match_rows.columns:
            r = match_rows.iloc[0]
            positions_plot.append(pos_aa)
            energies_plot.append(r["total energy"])

    if len(positions_plot) > 0:
        fig, ax = plt.subplots(figsize=(8, 3.5))
        ax.scatter(positions_plot, energies_plot, marker="x")
        ax.set_xlabel("Residue position")
        ax.set_ylabel("total energy")
        ax.set_title("FoldX total energy across sequence positions")
        ax.grid(True, linestyle="--", alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(OUT_DIR, "FoldX_total_energy_vs_position.png"), dpi=300, bbox_inches="tight")
        plt.close()

    foldx_idx = []
    for v in aa_variants:
        m = re.match(r"^([A-Z])(\d+)([A-Z])$", v)
        if not m:
            continue
        wt, pos_aa, mut = m.group(1), int(m.group(2)), m.group(3)
        pos_fx = pos_aa - 1
        match_rows = foldx[
            (foldx["position"] == pos_fx) &
            (foldx["wt_aa"] == wt) &
            (foldx["mut_aa"] == mut)
        ]
        if len(match_rows) == 1:
            foldx_idx.append(match_rows.index[0])

    foldx_154 = foldx.loc[foldx_idx].copy()

    foldx_corr_cols = [
        "SD",
        "total energy",
        "Backbone Hbond",
        "Sidechain Hbond",
        "Van der Waals",
        "Electrostatics",
        "Solvation Polar",
        "Solvation Hydrophobic",
        "Van der Waals clashes",
        "entropy sidechain",
        "entropy mainchain",
        "torsional clash",
        "backbone clash",
        "helix dipole",
        "electrostatic kon",
        "energy Ionisation",
    ]
    foldx_corr_cols = [c for c in foldx_corr_cols if c in foldx_154.columns]

    if len(foldx_corr_cols) > 1:
        corr_foldx = foldx_154[foldx_corr_cols].corr()
        fig, ax = plt.subplots(figsize=(8, 6))
        im = ax.imshow(corr_foldx.values, vmin=-1, vmax=1)
        ax.set_xticks(range(len(foldx_corr_cols)))
        ax.set_yticks(range(len(foldx_corr_cols)))
        ax.set_xticklabels(foldx_corr_cols, rotation=90)
        ax.set_yticklabels(foldx_corr_cols)
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label("Correlation coefficient")
        ax.set_title("Correlation Between FoldX Energy Terms")
        plt.tight_layout()
        plt.savefig(os.path.join(OUT_DIR, "FoldX_energy_correlation_154.png"), dpi=300, bbox_inches="tight")
        plt.close()

    corr_cols = [
        "ddg_prediction",
        "delta_stability_encom",
        "delta_vibrational_entropy",
        "mcsm",
        "sdm",
        "duet",
        "total energy",
        "Van der Waals",
        "Solvation Polar",
        "Solvation Hydrophobic",
        "Van der Waals clashes",
        "entropy sidechain",
        "entropy mainchain",
        "revel_revel",
        "phred",
    ]
    corr_cols = [c for c in corr_cols if c in final.columns]

    if len(corr_cols) > 1:
        corr_data = final[corr_cols].apply(pd.to_numeric, errors="coerce")
        corr_mat = corr_data.corr()

        fig, ax = plt.subplots(figsize=(8, 7))
        im2 = ax.imshow(corr_mat.values, vmin=-1, vmax=1)
        ax.set_xticks(range(len(corr_cols)))
        ax.set_yticks(range(len(corr_cols)))
        ax.set_xticklabels(corr_cols, rotation=90)
        ax.set_yticklabels(corr_cols)
        cbar2 = plt.colorbar(im2, ax=ax, fraction=0.046, pad=0.04)
        cbar2.set_label("Pearson correlation coefficient")
        ax.set_title("Correlation Between Structural Predictors, REVEL and CADD")
        plt.tight_layout()
        plt.savefig(os.path.join(OUT_DIR, "integrated_corr_heatmap.png"), dpi=300, bbox_inches="tight")
        plt.close()

    if "label_num" in final.columns and "revel_revel" in final.columns:
        labels = pd.to_numeric(final["label_num"], errors="coerce")
        label_groups = sorted(labels.dropna().unique())

        box_revel = [final.loc[labels == g, "revel_revel"].dropna() for g in label_groups]
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.boxplot(box_revel, tick_labels=label_groups, showfliers=True)
        ax.set_xlabel("label_num")
        ax.set_ylabel("REVEL score")
        ax.set_title("REVEL Score Distribution Across Clinical Labels")
        ax.yaxis.grid(True, linestyle="--", alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(OUT_DIR, "revel_by_label.png"), dpi=300, bbox_inches="tight")
        plt.close()

    if "label_num" in final.columns and "phred" in final.columns:
        labels = pd.to_numeric(final["label_num"], errors="coerce")
        label_groups = sorted(labels.dropna().unique())

        box_phred = [final.loc[labels == g, "phred"].dropna() for g in label_groups]
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.boxplot(box_phred, tick_labels=label_groups, showfliers=True)
        ax.set_xlabel("label_num")
        ax.set_ylabel("CADD phred")
        ax.set_title("CADD Phred Distribution Across Clinical Labels")
        ax.yaxis.grid(True, linestyle="--", alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(OUT_DIR, "cadd_by_label.png"), dpi=300, bbox_inches="tight")
        plt.close()

    if "label_num" in final.columns and {"revel_revel", "ddg_prediction"} <= set(final.columns):
        labels = pd.to_numeric(final["label_num"], errors="coerce")
        label_groups = sorted(labels.dropna().unique())

        fig, ax = plt.subplots(figsize=(6, 4))
        for g in label_groups:
            mask = labels == g
            ax.scatter(final.loc[mask, "revel_revel"], final.loc[mask, "ddg_prediction"], label=str(g), alpha=0.8)
        ax.axhline(0, linestyle="--", linewidth=1)
        ax.set_xlabel("REVEL score")
        ax.set_ylabel("ddg_prediction")
        ax.set_title("ddg_prediction vs REVEL score")
        ax.legend(title="label_num")
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(OUT_DIR, "ddg_vs_revel_scatter.png"), dpi=300, bbox_inches="tight")
        plt.close()

        fig, ax = plt.subplots(figsize=(6, 4))
        for g in label_groups:
            mask = labels == g
            ax.scatter(final.loc[mask, "revel_revel"], final.loc[mask, "phred"], label=str(g), alpha=0.8)
        ax.set_xlabel("REVEL score")
        ax.set_ylabel("CADD phred")
        ax.set_title("REVEL vs CADD phred by clinical label")
        ax.legend(title="label_num")
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(OUT_DIR, "revel_vs_phred_scatter.png"), dpi=300, bbox_inches="tight")
        plt.close()

    print("Integrated analysis finished.")


if __name__ == "__main__":
    main()