import pandas as pd
import numpy as np
from pathlib import Path
import configparser
import json
import sys
from typing import Sequence
import plotly.figure_factory as ff
import plotly
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from tqdm import tqdm

import statsmodels.api as sm
from scipy.stats import nbinom

def negbinom_fit(column:pd.Series):
    nb_fit = sm.NegativeBinomial(column, np.ones_like(column)).fit(start_params=[1,1], disp=0) #disp=0 == quiet 
    mu = np.exp(nb_fit.params[0])
    over_disp = nb_fit.params[1]
    p = 1/(1+mu*over_disp)
    n = mu*p/(1-p)
    return mu, over_disp, p, n

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--dir-count-matrices", required=True, type=Path)
    parser.add_argument("--dir-mag-screened", required=True, type=Path)
    parser.add_argument("--dataset", required=True, type=str)
    parser.add_argument("-o", required=True, type=Path)
    args = parser.parse_args()

    dir_count_matrices = Path(args.dir_count_matrices)   # Path("../data/simulated_data_large/kmer_quantification/count_matrices/")
    dir_mag_gene_sets = Path(args.dir_mag_screened)      #Path("../data/simulated_data_large/MAGpy/screened_flat/")

    dataset=args.dataset    #"0_5GB"

    outdir = Path(args.o) / dataset   #Path("../data/simulated_data_large/results/10_maginator_nb_fit") / dataset
    outdir.mkdir(parents=True, exist_ok=True)

    catalogue_count_files = [file for file in dir_count_matrices.glob("*.tsv") if not file.stem=="counts_all"]

    # Collect counts:
    catalouge_count_files = [file for file in dir_count_matrices.glob("*.tsv") if not file.stem=="counts_all"]
    df_count = pd.concat(
        pd.read_csv(file, sep="\t", index_col=0)\
            .reset_index()\
            .rename(columns={'index':'kmer'})\
            .melt(id_vars=["kmer"], var_name='dataset_sample', value_name='count')\
            .assign(catalogue_name = file.stem)
        for file in catalogue_count_files
    ).reset_index(drop=True)
    df_count[["dataset","sample"]] = df_count['dataset_sample'].str.split(".", expand=True)

    # Collect id sets
    id_frames = []
    for file in dir_mag_gene_sets.glob("*_kmers.csv"):
        cols = {"catalogue_name", "random","init","best"}
        df_i = pd.read_csv(file)\
            .assign(catalogue_name = file.stem.rsplit("_kmers",1)[0])
        cols_out = list(cols.intersection(set(df_i.columns)))
        id_frames.append(df_i[cols_out])
    df_id_sets = pd.concat(id_frames)

    # df_id_sets = pd.concat(
    #     pd.read_csv(file)\
    #         .assign(catalogue_name = file.stem.rsplit("_kmers",1)[0])[["catalogue_name", "random","init","best"]]
    #     for file in dir_mag_gene_sets.glob("*_kmers.csv")
    # )
    df_ids_long = df_id_sets.melt(
        id_vars=["catalogue_name"],
        value_vars = ["random","init","best"], 
        value_name = "kmer",
        var_name   = "catalogue_type"
    )

    df_gene_long_total = pd.concat(
        [
            df_ids_long,
            df_count[["kmer", 'catalogue_name']].drop_duplicates().assign(catalogue_type = "all_10000")
        ]
    ).reset_index(drop=True)

    #Generating figures:
    # We only look at a single dataset to reduce clutter
    df_count_dataset = df_count.query(f"dataset == '{dataset}'") 
     
    catalouges = df_count_dataset.catalogue_name.drop_duplicates()
    colormap = { #Yep this is wierd setup
        method: ["blue", "red","goldenrod", "green"][i]
        for i, method in enumerate(["all_10000","random","init","best"])
        }

    for cat in (pbar := tqdm(catalouges)):
        pbar.set_description(f"[{cat}] Wrangling data")
        df_cat_count = pd.merge(
            df_gene_long_total.query(f"catalogue_name == '{cat}'"),
            df_count_dataset.loc[:,["kmer","count"]],
            how="left",
            on = ["kmer"]
        )
        names  = []
        counts = []
        size   = []
        models = []
        for name, df_group in df_cat_count.groupby('catalogue_type'):
            if all(df_group.kmer.isna()):
                print(f"Skipping: [{cat}:{name}] - not defined in catalogue")
                continue

            names.append(name)
            counts_i = df_group["count"].tolist()
            counts.append(counts_i)
            size.append(len(counts_i))

            nb_fit = sm.NegativeBinomial(counts_i, np.ones_like(counts_i)).fit(disp=0, start_params=[1,1])
            models.append(nb_fit)
        
        pbar.set_description(f"[{cat}] Generating figure")
        # Create figure
        fig = make_subplots(
            rows=2, cols=1,
            shared_xaxes=True,
            row_width=[0.3, 0.7],
            #subplot_titles=[],
            vertical_spacing=0.06
        )
        
        colors = plotly.colors.qualitative.Safe
        for i in range(len(names)):
            fig.add_trace(
                go.Histogram(
                    legendgrouptitle_text = names[i],
                    name="Counts",
                    legendgroup = names[i],
                    x=counts[i],
                    histnorm='probability',
                    xbins = dict(size=1),
                    opacity=0.5,
                    marker = dict(color=colormap[names[i]])
                ),
            row=1, col=1
            )
            # Grabbing parameters from the fitted model
            mu = np.exp(models[i].params[0])
            over_disp = models[i].params[1]
            p = 1/(1+mu*over_disp)
            n = mu*p/(1-p)

            log_likelihood = models[i].llf
            max_count = int(max(counts[i]))
            x_fit = np.linspace(0, max_count, max_count+1)
            y_fit = nbinom.pmf(x_fit, n, p)

            fig.add_trace(
                go.Scatter(
                    x=x_fit,
                    y=y_fit,
                    legendgrouptitle_text = names[i],
                    mode="lines",
                    name="Fitted NB curve",
                    legendgroup = names[i],
                    marker = dict(color=colormap[names[i]])
                ),
            row=1, col=1
            )  
            
            # Adding rug to aid viewing the distributions.
            fig.add_trace(
                go.Scatter(
                    x = counts[i],
                    y = [names[i] for k in range(len(counts[i]))],
                    legendgroup = names[i],
                    name='rug',
                    mode='markers',
                    showlegend=False,
                    marker=dict(symbol='line-ns-open', color=colormap[names[i]], opacity=0.2)
                ),
                row=2, col=1
            )
            
        sub_title = "Fit LogLikelihoods: " + "; ".join([f"{name}: {int(model.llf)}" for name, model in zip(names, models)])
        sub_title2 = "Fit R2: " + "; ".join([f"{name}: {model.prsquared:.1e}" for name, model in zip(names, models)])
        fig.update_layout(title=f"Catalogue: [{cat}]<br>Normalised count historgram with fitted NB curve.") #<br><sup>{sub_title}<br>{sub_title2}")
        fig.update_layout(barmode="overlay")
        
        if cat == "NZ_CP053893.1.region002":
            x_max = 400
        else:
            shared_x_max = max(max(x) for x in counts)
            x_max = min(150, shared_x_max)
        
        fig.update_xaxes(range=[0,x_max])
        pbar.set_description(f"[{cat}] Writing figure")

        fig.write_image(outdir/(cat+".png"), width=800, height=650, scale=2)

        best_index = names.index("best")
        x_max_red = max(max(counts[best_index]), 10)
        fig.update_xaxes(range=[0, x_max_red])
        fig.write_image(outdir/ (cat+"_zoomed.png"))

