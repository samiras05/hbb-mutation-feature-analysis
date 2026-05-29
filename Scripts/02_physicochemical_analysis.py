import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "Data"

FILE_PATH = DATA_DIR / "HBB_10_merged_validated_ref_alt.xlsx"

FIG_DIR = "figures_physicochemical_analysis"
os.makedirs(FIG_DIR, exist_ok=True)

plt.rcParams["figure.dpi"] = 120
plt.rcParams["axes.titlesize"] = 11
plt.rcParams["axes.labelsize"] = 9
plt.rcParams["xtick.labelsize"] = 8
plt.rcParams["ytick.labelsize"] = 8


def main() -> None:
    df = pd.read_excel(FILE_PATH, header=1)

    numeric_cols = [
        "volume_ref", "volume_mut",
        "mass_ref", "mass_mut",
        "hydropathy_index_ref", "hydropathy_index_mut",
        "charge_ref", "charge_mut",
        "vol_delta", "charge_delta",
        "kd_wt", "kd_mt", "kd_delta",
        "mindist_to_fe_ang",
        "asa", "rsa"
    ]

    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    print("Number of variants:", df.shape[0])

    if {"volume_ref", "volume_mut", "mass_ref", "mass_mut"} <= set(df.columns):
        plt.figure(figsize=(10, 4))

        plt.subplot(1, 2, 1)
        plt.hist(df["volume_ref"].dropna(), bins=15, alpha=0.6, label="ref")
        plt.hist(df["volume_mut"].dropna(), bins=15, alpha=0.6, label="mut")
        plt.xlabel("Side-chain volume")
        plt.ylabel("Count")
        plt.title("Volume (ref vs mut)")
        plt.legend()

        plt.subplot(1, 2, 2)
        plt.hist(df["mass_ref"].dropna(), bins=15, alpha=0.6, label="ref")
        plt.hist(df["mass_mut"].dropna(), bins=15, alpha=0.6, label="mut")
        plt.xlabel("Side-chain mass")
        plt.ylabel("Count")
        plt.title("Mass (ref vs mut)")
        plt.legend()

        plt.tight_layout()
        plt.savefig(os.path.join(FIG_DIR, "volume_mass_ref_vs_mut.png"), dpi=300, bbox_inches="tight")
        plt.close()

    if {"hydropathy_index_ref", "hydropathy_index_mut"} <= set(df.columns):
        plt.figure(figsize=(5, 5))
        plt.scatter(df["hydropathy_index_ref"], df["hydropathy_index_mut"], alpha=0.7)
        plt.xlabel("Hydropathy (ref)")
        plt.ylabel("Hydropathy (mut)")
        plt.title("Hydropathy shift due to mutation")

        min_val = np.nanmin([df["hydropathy_index_ref"].min(), df["hydropathy_index_mut"].min()])
        max_val = np.nanmax([df["hydropathy_index_ref"].max(), df["hydropathy_index_mut"].max()])
        plt.plot([min_val, max_val], [min_val, max_val], linestyle="--")

        plt.tight_layout()
        plt.savefig(os.path.join(FIG_DIR, "hydropathy_shift_scatter.png"), dpi=300, bbox_inches="tight")
        plt.close()

    plt.figure(figsize=(10, 4))

    if "vol_delta" in df.columns:
        plt.subplot(1, 2, 1)
        plt.hist(df["vol_delta"].dropna(), bins=15)
        plt.xlabel("Δvolume")
        plt.ylabel("Count")
        plt.title("Distribution of Δvolume")

    if "charge_delta" in df.columns:
        plt.subplot(1, 2, 2)
        plt.hist(df["charge_delta"].dropna(), bins=9)
        plt.xlabel("Δcharge")
        plt.ylabel("Count")
        plt.title("Distribution of Δcharge")

    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, "delta_volume_charge_hist.png"), dpi=300, bbox_inches="tight")
    plt.close()

    if "bin_10aa" in df.columns:
        plt.figure(figsize=(8, 4))
        df["bin_10aa"].value_counts().sort_index().plot(kind="bar")
        plt.xlabel("10-aa bin")
        plt.ylabel("Number of variants")
        plt.title("Distribution of variants along the sequence")
        plt.tight_layout()
        plt.savefig(os.path.join(FIG_DIR, "variant_distribution_by_bin.png"), dpi=300, bbox_inches="tight")
        plt.close()

    if "ss" in df.columns:
        plt.figure(figsize=(6, 4))
        df["ss"].value_counts().plot(kind="bar")
        plt.xlabel("Secondary structure type")
        plt.ylabel("Number of variants")
        plt.title("Variants by Secondary Structure")

        annotation_text = (
            "H: alpha-Helix\n"
            "G: 3-10 Helix\n"
            "I: pi-Helix\n"
            "T: Turn\n"
            "S: Bend\n"
            "- : Coil/Loop"
        )

        plt.text(
            0.98, 0.98, annotation_text,
            transform=plt.gca().transAxes,
            ha="right",
            va="top",
            fontsize=9,
            bbox=dict(facecolor="white", alpha=0.65, edgecolor="gray")
        )

        plt.tight_layout()
        plt.savefig(os.path.join(FIG_DIR, "secondary_structure_distribution.png"), dpi=300, bbox_inches="tight")
        plt.close()

    if "rsa" in df.columns:
        plt.figure(figsize=(6, 4))
        plt.hist(df["rsa"].dropna(), bins=15)
        plt.xlabel("Relative solvent accessibility (RSA)")
        plt.ylabel("Count")
        plt.title("Distribution of RSA for variant positions")
        plt.tight_layout()
        plt.savefig(os.path.join(FIG_DIR, "rsa_distribution.png"), dpi=300, bbox_inches="tight")
        plt.close()

    if {"mindist_to_fe_ang", "kd_delta"} <= set(df.columns):
        plt.figure(figsize=(6, 4))
        plt.scatter(df["mindist_to_fe_ang"], df["kd_delta"], alpha=0.7)
        plt.xlabel("Minimum distance to Fe (Å)")
        plt.ylabel("ΔKd")
        plt.title("Distance to heme vs ΔKd")
        plt.tight_layout()
        plt.savefig(os.path.join(FIG_DIR, "heme_distance_vs_kd_delta.png"), dpi=300, bbox_inches="tight")
        plt.close()

    if {"kd_wt", "kd_mt", "kd_delta"} <= set(df.columns):
        plt.figure(figsize=(10, 4))

        plt.subplot(1, 2, 1)
        plt.hist(df["kd_wt"].dropna(), bins=15, alpha=0.7, label="kd_wt")
        plt.hist(df["kd_mt"].dropna(), bins=15, alpha=0.7, label="kd_mt")
        plt.xlabel("Kd")
        plt.ylabel("Count")
        plt.title("Kd (wild-type vs mutant)")
        plt.legend()

        plt.subplot(1, 2, 2)
        plt.hist(df["kd_delta"].dropna(), bins=15)
        plt.xlabel("ΔKd")
        plt.ylabel("Count")
        plt.title("Distribution of ΔKd")

        plt.tight_layout()
        plt.savefig(os.path.join(FIG_DIR, "kd_distributions.png"), dpi=300, bbox_inches="tight")
        plt.close()

    print("Physicochemical analysis finished.")


if __name__ == "__main__":
    main()