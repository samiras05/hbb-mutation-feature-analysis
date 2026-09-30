import os
import re
from collections import Counter

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "Data"

FILE_PATH = DATA_DIR / "HBB_10_merged_validated_ref_alt.xlsx"

OUT_DIR = "dynamut_figures"
os.makedirs(OUT_DIR, exist_ok=True)

plt.rcParams["figure.dpi"] = 120
plt.rcParams["axes.titlesize"] = 11
plt.rcParams["axes.labelsize"] = 10
plt.rcParams["xtick.labelsize"] = 8
plt.rcParams["ytick.labelsize"] = 8




def main() -> None:
    df = pd.read_excel(FILE_PATH, header=1)

    dyn_cols = [
        "delta_vibrational_entropy",
        "ddg_prediction",
        "delta_stability_encom",
        "mcsm",
        "sdm",
        "duet",
    ]
    meta_cols = ["variant", "label_text", "label_num", "mutation"]

    present = [c for c in dyn_cols + meta_cols if c in df.columns]
    df_dyn = df[present].copy()

    for col in dyn_cols:
        if col in df_dyn.columns:
            df_dyn[col] = pd.to_numeric(df_dyn[col], errors="coerce")

    if "label_num" in df_dyn.columns:
        df_dyn["label_num"] = pd.to_numeric(df_dyn["label_num"], errors="coerce")

    if "variant" in df_dyn.columns:
        def extract_position(variant):
            m = re.search(r"(\d+)", str(variant))
            return int(m.group(1)) if m else np.nan

        df_dyn["position"] = df_dyn["variant"].apply(extract_position)

    print("=== Descriptive statistics for DynaMut columns ===")
    print(df_dyn[[c for c in dyn_cols if c in df_dyn.columns]].describe())

    for col in [c for c in dyn_cols if c in df_dyn.columns]:
        values = df_dyn[col].dropna()
        plt.figure(figsize=(6, 4))
        plt.hist(values, bins=18)
        plt.title(f"Distribution of {col}")
        plt.xlabel(col)
        plt.ylabel("Count")
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(OUT_DIR, f"dyn_hist_{col}.png"), bbox_inches="tight")
        plt.close()

    corr_matrix = df_dyn[[c for c in dyn_cols if c in df_dyn.columns]].corr(method="pearson")
    print("\n=== Correlation matrix between DynaMut predictors ===")
    print(corr_matrix)

    plt.figure(figsize=(6, 5))
    im = plt.imshow(corr_matrix.values)
    plt.xticks(range(len(corr_matrix.columns)), corr_matrix.columns, rotation=45, ha="right")
    plt.yticks(range(len(corr_matrix.index)), corr_matrix.index)
    plt.colorbar(im, fraction=0.046, pad=0.04)
    plt.title("Correlation matrix of DynaMut outputs")
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "dyn_corr_matrix.png"), bbox_inches="tight")
    plt.close()

    if {"ddg_prediction", "mcsm"} <= set(df_dyn.columns):
        pair_pairs = [
            ("ddg_prediction", "mcsm"),
            ("ddg_prediction", "duet"),
            ("mcsm", "sdm"),
        ]

        for x_col, y_col in pair_pairs:
            if x_col in df_dyn.columns and y_col in df_dyn.columns:
                plt.figure(figsize=(6, 4))
                plt.scatter(df_dyn[x_col], df_dyn[y_col], s=25, alpha=0.75)
                plt.xlabel(x_col)
                plt.ylabel(y_col)
                plt.title(f"{y_col} vs {x_col}")
                plt.grid(True, alpha=0.3)
                plt.tight_layout()
                plt.savefig(os.path.join(OUT_DIR, f"dyn_scatter_{y_col}_vs_{x_col}.png"), bbox_inches="tight")
                plt.close()

    sign_cols = [c for c in ["ddg_prediction", "delta_stability_encom", "mcsm", "sdm", "duet"] if c in df_dyn.columns]
    if len(sign_cols) >= 3:
        sign_df = np.sign(df_dyn[sign_cols])
        n_pos = (sign_df > 0).sum(axis=1)
        n_neg = (sign_df < 0).sum(axis=1)

        consensus_categories = []
        for pos_count, neg_count in zip(n_pos, n_neg):
            if neg_count == len(sign_cols):
                consensus_categories.append("all_destabilizing")
            elif pos_count == len(sign_cols):
                consensus_categories.append("all_stabilizing")
            elif neg_count >= max(4, len(sign_cols) - 1) and pos_count == 0:
                consensus_categories.append("mostly_destabilizing")
            elif pos_count >= max(4, len(sign_cols) - 1) and neg_count == 0:
                consensus_categories.append("mostly_stabilizing")
            else:
                consensus_categories.append("mixed")

        df_dyn["consensus_category"] = consensus_categories
        consensus_counts = Counter(consensus_categories)

        print("\n=== Consensus categories ===")
        print(consensus_counts)

        plt.figure(figsize=(6, 4))
        plt.bar(list(consensus_counts.keys()), list(consensus_counts.values()))
        plt.xlabel("Consensus category")
        plt.ylabel("Number of variants")
        plt.title("Consensus of stability prediction among DynaMut tools")
        plt.xticks(rotation=20, ha="right")
        plt.grid(True, axis="y", alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(OUT_DIR, "dyn_consensus_categories.png"), bbox_inches="tight")
        plt.close()

    if "delta_vibrational_entropy" in df_dyn.columns and "ddg_prediction" in df_dyn.columns:
        entropy_col = "delta_vibrational_entropy"
        entropy_corr = df_dyn[[entropy_col] + [c for c in dyn_cols if c in df_dyn.columns and c != entropy_col]].corr().loc[entropy_col]
        print("\n=== Correlation of delta_vibrational_entropy with other DynaMut outputs ===")
        print(entropy_corr)

        plt.figure(figsize=(6, 4))
        plt.scatter(df_dyn[entropy_col], df_dyn["ddg_prediction"], s=25, alpha=0.75)
        plt.axhline(0, linewidth=1)
        plt.axvline(0, linewidth=1)
        plt.xlabel("delta_vibrational_entropy")
        plt.ylabel("ddg_prediction")
        plt.title("ddg_prediction vs delta_vibrational_entropy")
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(OUT_DIR, "dyn_scatter_ddg_vs_entropy.png"), bbox_inches="tight")
        plt.close()

    if "position" in df_dyn.columns and "ddg_prediction" in df_dyn.columns:
        pos_stats = (
            df_dyn.groupby("position")["ddg_prediction"]
            .agg(["count", "mean", "min", "max"])
            .dropna()
            .sort_index()
        )

        print("\n=== Position-based ddg_prediction summary (head) ===")
        print(pos_stats.head(10))

        plt.figure(figsize=(7, 4))
        plt.plot(pos_stats.index, pos_stats["mean"], linewidth=1.5)
        plt.xlabel("Residue position")
        plt.ylabel("Mean ddg_prediction")
        plt.title("Mean ddg_prediction along HBB sequence")
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(OUT_DIR, "dyn_position_mean_ddg.png"), bbox_inches="tight")
        plt.close()

    if "ddg_prediction" in df_dyn.columns:
        sorted_df = df_dyn.sort_values("ddg_prediction")
        print("\n=== Top 10 most destabilizing variants ===")
        print(sorted_df.head(10)[[c for c in ["variant", "ddg_prediction", "label_num", "label_text"] if c in sorted_df.columns]])

        print("\n=== Top 10 most stabilizing variants ===")
        print(sorted_df.tail(10)[[c for c in ["variant", "ddg_prediction", "label_num", "label_text"] if c in sorted_df.columns]])

        top_destab = sorted_df.head(10)
        plt.figure(figsize=(7, 4))
        plt.bar(top_destab["variant"].astype(str), top_destab["ddg_prediction"])
        plt.xlabel("Variant")
        plt.ylabel("ddg_prediction")
        plt.title("Top 10 most destabilizing variants")
        plt.xticks(rotation=45, ha="right")
        plt.grid(True, axis="y", alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(OUT_DIR, "dyn_top10_destabilizing.png"), bbox_inches="tight")
        plt.close()

    if all(c in df_dyn.columns for c in ["ddg_prediction", "delta_stability_encom", "mcsm", "sdm", "duet"]):
        pca_cols = ["ddg_prediction", "delta_stability_encom", "mcsm", "sdm", "duet"]
        tmp = df_dyn[pca_cols].dropna().copy()
        if len(tmp) > 2:
            scaler = StandardScaler()
            X = scaler.fit_transform(tmp)

            pca = PCA(n_components=2, random_state=42)
            pcs = pca.fit_transform(X)

            plt.figure(figsize=(6, 5))
            if "label_num" in df_dyn.columns:
                labels = df_dyn.loc[tmp.index, "label_num"]
                scatter = plt.scatter(pcs[:, 0], pcs[:, 1], c=labels, cmap="viridis", alpha=0.8)
                plt.colorbar(scatter, label="label_num")
            else:
                plt.scatter(pcs[:, 0], pcs[:, 1], alpha=0.8)

            plt.xlabel("PC1")
            plt.ylabel("PC2")
            plt.title("PCA of DynaMut predictors")
            plt.tight_layout()
            plt.savefig(os.path.join(OUT_DIR, "dyn_pca.png"), bbox_inches="tight")
            plt.close()

            if len(tmp) >= 3:
                km = KMeans(n_clusters=min(3, len(tmp)), random_state=42, n_init=10)
                clusters = km.fit_predict(X)

                plt.figure(figsize=(6, 5))
                plt.scatter(pcs[:, 0], pcs[:, 1], c=clusters, cmap="tab10", alpha=0.8)
                plt.xlabel("PC1")
                plt.ylabel("PC2")
                plt.title("KMeans clustering on DynaMut predictors")
                plt.tight_layout()
                plt.savefig(os.path.join(OUT_DIR, "dyn_kmeans_clusters.png"), bbox_inches="tight")
                plt.close()

    print("DynaMut analysis finished.")


if __name__ == "__main__":
    main()