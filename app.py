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


def plot_volcano(rd, comparison, annot_col="name", gene_list=None):
    df = pd.concat([rd.annotation, rd.comparisons[comparison].df], axis=1, join="inner")

    st.subheader("Volcano plot")
    plt.figure()

    v = Volcano(df, annot_col=annot_col)

    fig, axes = plt.subplots(1, 1, figsize=(8, 6))
    v.plot()

    if gene_list:
        sel = df.query("@annot_col in @gene_list").copy()

        sel.loc[:, "y"] = [20, 25, 120]
        sel.loc[:, "x"] = [2, 3.5, -3]

        plt.scatter(sel.log2FoldChange, -np.log10(sel.padj), facecolor="none", edgecolor="black")
        sel.apply(plot_gene_names, axis=1)

    plt.ylabel("-log10(adjusted pvalue)")
    plt.xlabel("Fold Change")
    plt.xlim([-6, 6])

    st.pyplot(plt)


def plot_heatmap(rd, comparison_sel, regulation, condition, groups, annot_col="name", log2FC=0, alpha=0.05):
    # gene_ids = rd.annotation.query("gene_name in @total_gene_list[@gene_list]").index

    # rd.log2_fc = log2FC
    # rd.alpha = alpha

    gene_list = rd.get_gene_lists()[comparison_sel][regulation]
    print(f"Gene list after removing genes not significative N: {len(gene_list)}")

    count_indexes = rd.counts_norm.index

    gene_list = set(gene_list).intersection(set(count_indexes))
    print(f"Gene list after removing genes not present in counts matrix N: {len(gene_list)}")

    samples = list(rd.design_df.query(f"{condition} in @groups").index)

    data = rd.counts_norm.loc[list(gene_list), samples]
    # data = rd.counts_norm.loc[list(gene_list),:]

    annotation = rd.df.loc[:, ("annotation", annot_col)].rename(index={"f(annotation, {annot_col})": annot_col})
    annotation.name = annot_col
    data = pd.concat([data, annotation], axis=1, join="inner").set_index(annot_col)

    # Removing genes with expression standard deviation == 0 (otherwise z_score cannot be computed)
    print(f"Number of genes with null std removed: {(data.std(axis=1) == 0).sum()}")
    data = data.loc[data.std(axis=1) != 0, :]

    sns.set(font_scale=1)
    p = Clustermap(
        data,
        sample_groups_df=rd.design_df.query(f"{condition} in @groups"),
        z_score=0,
        center=0,
        sample_groups_sel=[condition],
    )
    p.params["legend.sample.bbox"] = (10, 2)
    p.plot()


def main():
    st.title("Interactive RNADiff")

    st.sidebar.subheader("Input file")
    rd_file = st.sidebar.file_uploader("Upload RNADiffResults pickled object:")

    if rd_file:
        rd = pickle.load(rd_file)

        st.subheader("Design")
        st.write(rd.design_df)

        st.subheader("Results")
        st.write(rd.df)

        st.sidebar.subheader("General parameters:")
        comparison = st.sidebar.selectbox("Pick comparison:", rd.get_gene_lists().keys())
        regulation = st.sidebar.selectbox("Pick regulation:", ["down", "up", "all", "no_subsetting"])
        annot_col = st.sidebar.selectbox("Pick annotation column:", rd.df.annotation.columns)

        st.sidebar.subheader("Heatmap options:")
        condition = st.sidebar.selectbox("Pick condition:", rd.design_df.columns.unique())
        groups = st.sidebar.multiselect("Pick groups:", rd.design_df[condition].unique())
        # log2fc = st.sidebar.number_input("Log2FC:", 1)
        # alpha = st.sidebar.number_input("padj threshold:", 0.05)

        plot_volcano(rd, comparison, annot_col=annot_col)

        plt.figure()
        plot_heatmap(rd, comparison, regulation, condition, groups, annot_col=annot_col)  # , log2fc, alpha)
        st.pyplot(plt)


if __name__ == "__main__":
    main()
