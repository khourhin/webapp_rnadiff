import streamlit as st
import pickle
from sequana.viz.volcano import Volcano
from sequana.viz.heatmap import Clustermap
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns

st.set_page_config(page_title="Heatmaps", page_icon="imgs/logo_256x256.png", layout="wide")


def plot_gene_names(row):
    plt.annotate(
        row.gene_name,
        (row.log2FoldChange, -np.log10(row.padj)),
        (row.log2FoldChange - 1, -np.log10(row.padj) - 1),
        arrowprops=dict(arrowstyle="->", color="black"),
    )


def plot_volcano(df, gene_list):
    v = Volcano(df, annot_col="gene_name")
    sel = df.query("gene_name in @gene_list").copy()

    sel.loc[:, "y"] = [20, 25, 120]
    sel.loc[:, "x"] = [2, 3.5, -3]

    fig, axes = plt.subplots(1, 1, figsize=(8, 6))
    v.plot()
    plt.scatter(sel.log2FoldChange, -np.log10(sel.padj), facecolor="none", edgecolor="black")
    sel.apply(plot_gene_names, axis=1)
    plt.ylabel("-log10(adjusted pvalue)")
    plt.xlabel("Fold Change")
    plt.xlim([-6, 6])


def volcano_plot(rd):
    df = pd.concat([rd.annotation, rd.comparisons["MTFP1_KO_vs_WT"].df], axis=1, join="inner")

    st.subheader("Volcano plot")
    plt.figure()
    plot_volcano(df, ["Nppa", "Nppb", "Mtfp1"])
    st.pyplot(plt)


def plot_heatmap(rd, comparison_sel, regulation, conditions, log2FC, alpha):
    # gene_ids = rd.annotation.query("gene_name in @total_gene_list[@gene_list]").index

    rd.log2_fc = log2FC
    rd.alpha = alpha

    gene_list = rd.get_gene_lists()[comparison_sel][regulation]
    print(f"Gene list after removing genes not significative N: {len(gene_list)}")

    count_indexes = rd.counts_norm.index

    gene_list = set(gene_list).intersection(set(count_indexes))
    print(f"Gene list after removing genes not present in counts matrix N: {len(gene_list)}")

    samples = list(rd.design_df.query("condition in @conditions").index)
    data = rd.counts_norm.loc[list(gene_list), samples]
    # data = rd.counts_norm.loc[list(gene_list),:]

    data = pd.concat([data, rd.annotation.gene_name], axis=1, join="inner").set_index("gene_name")

    # Removing genes with expression standard deviation == 0 (otherwise z_score cannot be computed)
    print(f"Number of genes with null std removed: {(data.std(axis=1) == 0).sum()}")
    data = data.loc[data.std(axis=1) != 0, :]

    sns.set(font_scale=1)
    p = Clustermap(
        data,
        sample_groups_df=rd.design_df.query("condition in @conditions"),
        z_score=0,
        center=0,
        sample_groups_sel=["condition"],
    )
    p.params["legend.sample.bbox"] = (10, 2)
    p.plot()


def main():
    st.title("Heatmaps")

    rd_file = st.sidebar.file_uploader("Upload RMADiffResults pickled object:")

    if rd_file:
        rd = pickle.load(rd_file)

        st.subheader("Design")
        st.write(rd.design_df)

        st.subheader("Results")
        st.write(rd.df)

        volcano_plot(rd)

        comparison = st.sidebar.selectbox("Pick comparison:", rd.get_gene_lists().keys())
        regulation = st.sidebar.selectbox("Pick regulation:", ["down", "up", "all", "no_subsetting"])
        conditions = st.sidebar.multiselect("Pick conditions:", rd.design_df["condition"].unique())
        log2fc = st.sidebar.number_input("Log2FC:", 1)
        alpha = st.sidebar.number_input("padj threshold:", 0.05)

        plt.figure()
        plot_heatmap(rd, comparison, regulation, conditions, log2fc, alpha)
        st.pyplot(plt)

        comparison_list = st.sidebar.multiselect("Pick comparisons:", rd.get_gene_lists().keys())


if __name__ == "__main__":
    main()
