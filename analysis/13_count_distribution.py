import sys

from pathlib import Path
import pandas as pd


from scipy.stats import nbinom
import statsmodels.api as sm
import numpy as np

import plotly.express as px
import plotly.graph_objects as go
import plotly.figure_factory as pff

import plotly.io as pio
scope = pio.kaleido.scope

scope._shutdown_kaleido() #attempt to fix kaleido hanging

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--dir-count-matrices", type=Path, required=True)
parser.add_argument("-o", type=Path, required=True)
args = parser.parse_args()


counts_dir = args.dir_count_matrices    #Path("data/simulated_data_large/kmer_quantification/count_matrices/")
glob_str = "*.tsv"
cluster_counts = counts_dir.glob(glob_str)
cluster_list = list(x for x in cluster_counts if not x.name.startswith("counts_all"))
outdir = args.o #Path("data/simulated_data_large/results/13_count_distribution")
outdir_by_dataset = outdir/"datasets"
outdir_by_catalogue = outdir/"catalogue"

outdir_by_dataset.mkdir(parents=True, exist_ok=True)
outdir_by_catalogue.mkdir(parents=True, exist_ok=True)


all_counts = []

if not cluster_list:
    raise FileNotFoundError(f"{counts_dir}/{glob_str}")

print("Generating per cluster figures", file=sys.stderr)

for fp_count in cluster_list:
    cluster_name = fp_count.stem
    df_count = pd.read_csv(fp_count, sep="\t", index_col=0)

    df_plot = df_count.reset_index().rename(columns={'index':'kmer'}).melt(id_vars="kmer", var_name="sample", value_name="counts")
    df_plot[["dataset","sample"]] = df_plot["sample"].str.split(".", expand=True)
    all_counts.append(df_plot.assign(cluster = cluster_name))

    values = []
    ds = []
    for name, df_ds in df_plot.groupby("dataset"):
        ds.append(name)
        values.append(df_ds["counts"].values.tolist())
    fig = pff.create_distplot(
        values, ds, bin_size=1
    )
    fig.update_layout(
        xaxis_title = "Counts",
        yaxis_title = "Probability Density",
        title=f"BGC_familiy: {cluster_name} per dataset. <br><sup>Probability density historgram with KDE curves.",
        yaxis2=dict(range=[0,2])
    )
    
    #fig.update_xaxes(range=[0, 250])
    outfile = (outdir_by_catalogue / cluster_name).with_suffix(".png")
    fig.write_image(outfile, height=1000, width=1000)
    scope._shutdown_kaleido()
    

print("Generating per dataset figures", file=sys.stderr)
df_total = pd.concat(all_counts)

for ds_name, df_ds in df_total.groupby("dataset"):
    
    values = []
    labels = []
    for label, df_i in df_ds.groupby("cluster"):
        labels.append(label)
        values.append(df_i["counts"].values.tolist())

    fig = pff.create_distplot(
        values, labels
    )
    fig.update_layout(
        xaxis_title = "Counts",
        yaxis_title = "Probability Density",
        title=f"Dataset: {ds_name} for all clusters.<br><sup>Probability density historgram with KDE curves."
    )
    fig.update_yaxes(range=[0,2])
    outfile = (outdir_by_dataset / ds_name).with_suffix(".png")
    fig.write_image(outfile, height=1000, width=1000)
    scope._shutdown_kaleido()
    

