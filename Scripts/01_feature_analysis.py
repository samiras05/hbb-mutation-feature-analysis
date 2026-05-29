import os
import pandas as pd
import matplotlib.pyplot as plt

from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "Data"

FILE_PATH = DATA_DIR / "HBB_10_merged_validated_ref_alt.xlsx"

FIG_DIR = "figures_feature_analysis"
os.makedirs(FIG_DIR, exist_ok=True)

plt.rcParams["figure.dpi"] = 120
plt.rcParams["axes.titlesize"] = 11
plt.rcParams["axes.labelsize"] = 9
plt.rcParams["xtick.labelsize"] = 8
plt.rcParams["ytick.labelsize"] = 8


def main() -> None:
    df = pd.read_excel(FILE_PATH, header=1)

    numeric_cols = ["pos_norm", "grantham", "blosum62", "label_num"]
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    needed = ["variant", "label_text", "label_num", "pos_norm", "wt_group", "mut_group", "is_same_group", "grantham", "blosum62"]
    present = [c for c in needed if c in df.columns]
    print("Head of key columns:")
    print(df[present].head())

    if "pos_norm" in df.columns:
        print("\n=== pos_norm descriptive statistics ===")
        print(df["pos_norm"].describe())

        plt.figure()
        df["pos_norm"].hist(bins=20)
        plt.xlabel("pos_norm")
        plt.ylabel("Count")
        plt.title("Distribution of pos_norm")
        plt.tight_layout()
        plt.savefig(os.path.join(FIG_DIR, "pos_norm_histogram.png"), dpi=300, bbox_inches="tight")
        plt.close()

        if "label_num" in df.columns:
            plt.figure()
            df.boxplot(column="pos_norm", by="label_num")
            plt.xlabel("label_num")
            plt.ylabel("pos_norm")
            plt.title("pos_norm by label_num")
            plt.suptitle("")
            plt.tight_layout()
            plt.savefig(os.path.join(FIG_DIR, "pos_norm_by_labelnum_boxplot.png"), dpi=300, bbox_inches="tight")
            plt.close()

    if "wt_group" in df.columns:
        print("\n=== wt_group counts ===")
        wt_counts = df["wt_group"].value_counts()
        print(wt_counts)

        plt.figure()
        wt_counts.plot(kind="bar")
        plt.xlabel("wt_group")
        plt.ylabel("Count")
        plt.title("Distribution of wt_group")
        plt.tight_layout()
        plt.savefig(os.path.join(FIG_DIR, "wt_group_barplot.png"), dpi=300, bbox_inches="tight")
        plt.close()

    if "mut_group" in df.columns:
        print("\n=== mut_group counts ===")
        mut_counts = df["mut_group"].value_counts()
        print(mut_counts)

        plt.figure()
        mut_counts.plot(kind="bar")
        plt.xlabel("mut_group")
        plt.ylabel("Count")
        plt.title("Distribution of mut_group")
        plt.tight_layout()
        plt.savefig(os.path.join(FIG_DIR, "mut_group_barplot.png"), dpi=300, bbox_inches="tight")
        plt.close()

    if "is_same_group" in df.columns:
        print("\n=== is_same_group counts ===")
        same_group_counts = df["is_same_group"].value_counts().sort_index()
        print(same_group_counts)

        plt.figure()
        same_group_counts.plot(kind="bar")
        plt.xlabel("is_same_group")
        plt.ylabel("Count")
        plt.title("Distribution of is_same_group")
        plt.tight_layout()
        plt.savefig(os.path.join(FIG_DIR, "is_same_group_barplot.png"), dpi=300, bbox_inches="tight")
        plt.close()

    if "wt_group" in df.columns and "mut_group" in df.columns:
        group_table = pd.crosstab(df["wt_group"], df["mut_group"])
        print("\n=== Crosstab wt_group vs mut_group ===")
        print(group_table)

        plt.figure()
        plt.imshow(group_table.values)
        plt.xticks(range(len(group_table.columns)), group_table.columns, rotation=45)
        plt.yticks(range(len(group_table.index)), group_table.index)
        plt.xlabel("mut_group")
        plt.ylabel("wt_group")
        plt.title("wt_group vs mut_group counts")
        plt.colorbar()
        plt.tight_layout()
        plt.savefig(os.path.join(FIG_DIR, "wt_vs_mut_group_heatmap.png"), dpi=300, bbox_inches="tight")
        plt.close()

    if "grantham" in df.columns:
        print("\n=== Grantham score descriptive statistics ===")
        print(df["grantham"].describe())

        plt.figure()
        df["grantham"].hist(bins=20)
        plt.xlabel("Grantham score")
        plt.ylabel("Count")
        plt.title("Distribution of Grantham scores")
        plt.tight_layout()
        plt.savefig(os.path.join(FIG_DIR, "grantham_histogram.png"), dpi=300, bbox_inches="tight")
        plt.close()

        if "label_num" in df.columns:
            plt.figure()
            df.boxplot(column="grantham", by="label_num")
            plt.xlabel("label_num")
            plt.ylabel("Grantham score")
            plt.title("Grantham score by label_num")
            plt.suptitle("")
            plt.tight_layout()
            plt.savefig(os.path.join(FIG_DIR, "grantham_by_labelnum_boxplot.png"), dpi=300, bbox_inches="tight")
            plt.close()

    if "blosum62" in df.columns:
        print("\n=== BLOSUM62 descriptive statistics ===")
        print(df["blosum62"].describe())

        plt.figure()
        df["blosum62"].hist(bins=15)
        plt.xlabel("BLOSUM62 score")
        plt.ylabel("Count")
        plt.title("Distribution of BLOSUM62 scores")
        plt.tight_layout()
        plt.savefig(os.path.join(FIG_DIR, "blosum62_histogram.png"), dpi=300, bbox_inches="tight")
        plt.close()

        if "label_num" in df.columns:
            plt.figure()
            df.boxplot(column="blosum62", by="label_num")
            plt.xlabel("label_num")
            plt.ylabel("BLOSUM62 score")
            plt.title("BLOSUM62 score by label_num")
            plt.suptitle("")
            plt.tight_layout()
            plt.savefig(os.path.join(FIG_DIR, "blosum62_by_labelnum_boxplot.png"), dpi=300, bbox_inches="tight")
            plt.close()

    if "grantham" in df.columns and "blosum62" in df.columns:
        plt.figure()
        plt.scatter(df["grantham"], df["blosum62"], alpha=0.7)
        plt.xlabel("Grantham score")
        plt.ylabel("BLOSUM62 score")
        plt.title("Grantham vs BLOSUM62")
        plt.tight_layout()
        plt.savefig(os.path.join(FIG_DIR, "grantham_vs_blosum62_scatter.png"), dpi=300, bbox_inches="tight")
        plt.close()

    print("\nAll analyses finished. Figures are saved in:", FIG_DIR)


if __name__ == "__main__":
    main()