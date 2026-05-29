import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "Data"

FILE_PATH = DATA_DIR / "FoldX.xlsx"

OUT_DIR = "foldx_figures"
os.makedirs(OUT_DIR, exist_ok=True)

plt.rcParams["figure.dpi"] = 120
plt.rcParams["axes.titlesize"] = 11
plt.rcParams["axes.labelsize"] = 9
plt.rcParams["xtick.labelsize"] = 8
plt.rcParams["ytick.labelsize"] = 8


def main() -> None:
    foldx = pd.read_excel(FILE_PATH)

    energy_cols_all = [
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
    "sloop_entropy",
    "mloop_entropy",
    "torsional clash",
    "backbone clash",
    "helix dipole",
    "water bridge",
    "electrostatic kon",
    "energy Ionisation",
]
    for col in energy_cols_all:
        foldx[col] = pd.to_numeric(foldx[col], errors="coerce")

    stds = foldx[energy_cols_all].describe().loc["std"]
    non_const_cols = stds[stds > 0].index.tolist()
    const_cols = stds[stds == 0].index.tolist()

    print("Non-constant energy terms:", non_const_cols)
    print("Constant (all ~0) terms:", const_cols)

    n = len(non_const_cols)
    if n > 0:
        n_cols = 4
        n_rows = int(np.ceil(n / n_cols))
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(12, 3 * n_rows))
        axes = np.array(axes).flatten()

        for ax, col in zip(axes, non_const_cols):
            vals = pd.to_numeric(foldx[col], errors="coerce").dropna()
            ax.hist(vals, bins=18)
            ax.set_title(col)
            ax.set_xlabel(col)
            ax.set_ylabel("Count")
            ax.grid(True, alpha=0.2)

        for ax in axes[n:]:
            ax.axis("off")

        fig.suptitle("Distributions of FoldX energy terms", fontsize=11)
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        plt.savefig(os.path.join(OUT_DIR, "FoldX_energy_histograms.png"), bbox_inches="tight")
        plt.close()

    if len(non_const_cols) > 1:
        corr_foldx = foldx[non_const_cols].corr()

        fig, ax = plt.subplots(figsize=(8, 6))
        im = ax.imshow(corr_foldx.values, vmin=-1, vmax=1)
        ax.set_xticks(range(len(non_const_cols)))
        ax.set_yticks(range(len(non_const_cols)))
        ax.set_xticklabels(non_const_cols, rotation=90)
        ax.set_yticklabels(non_const_cols)
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label("Correlation coefficient")
        ax.set_title("Correlation Between FoldX Energy Terms")
        plt.tight_layout()
        plt.savefig(os.path.join(OUT_DIR, "FoldX_energy_correlation.png"), dpi=300, bbox_inches="tight")
        plt.close()

        if "total energy" in foldx.columns:
            fig, ax = plt.subplots(figsize=(8, 3.5))
            x = np.arange(len(foldx))
            if "position" in foldx.columns:
                x = pd.to_numeric(foldx["position"], errors="coerce")
            ax.scatter(x, foldx["total energy"], marker="x")
            ax.set_xlabel("Residue position")
            ax.set_ylabel("total energy")
            ax.set_title("FoldX total energy across sequence positions")
            ax.grid(True, linestyle="--", alpha=0.3)
            plt.tight_layout()
            plt.savefig(os.path.join(OUT_DIR, "FoldX_total_energy_vs_position.png"), dpi=300, bbox_inches="tight")
            plt.close()

        tmp = foldx[non_const_cols].dropna()
        if len(tmp) > 2:
            X = StandardScaler().fit_transform(tmp)
            pcs = PCA(n_components=2, random_state=42).fit_transform(X)

            fig, ax = plt.subplots(figsize=(6, 5))
            ax.scatter(pcs[:, 0], pcs[:, 1], alpha=0.8)
            ax.set_xlabel("PC1")
            ax.set_ylabel("PC2")
            ax.set_title("PCA of FoldX energy terms")
            plt.tight_layout()
            plt.savefig(os.path.join(OUT_DIR, "FoldX_pca.png"), dpi=300, bbox_inches="tight")
            plt.close()

    print("FoldX analysis finished.")


if __name__ == "__main__":
    main()