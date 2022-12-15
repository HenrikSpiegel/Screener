import sys

from pathlib import Path
import pandas as pd


from scipy.stats import nbinom, poisson
import statsmodels.api as sm
import numpy as np

import plotly.graph_objects as go
from plotly.subplots import make_subplots

import plotly.io as pio
scope = pio.kaleido.scope

scope._shutdown_kaleido() #attempt to fix kaleido hanging

def is_outlier_iqr(values: pd.Series, iqr_weight=1.5):
    q1, q3 = values.quantile([0.25, 0.75]).values.flatten()
    IQR = q3-q1
    outlier_dist = iqr_weight*IQR
    
    return (values < q1-outlier_dist) | (values > q3+outlier_dist)
    
    

def negbinom_fit(column:pd.Series):
    nb_fit = sm.NegativeBinomial(column, np.ones_like(column)).fit(start_params=[1,1], disp=0) #disp=0 == quiet 
    mu = np.exp(nb_fit.params[0])
    over_disp = nb_fit.params[1]
    p = 1/(1+mu*over_disp)
    n = mu*p/(1-p)
    return mu, over_disp, p, n


# def poisson_fit(counts) -> float:
#     pois_fit = sm.Poisson(counts, np.ones_like(counts)).fit(start_params=[1], disp=0)
#     mu = np.exp(pois_fit.params[0])
#     return mu


def generate_NB_trace(counts):
    mu, over_disp, p, n = negbinom_fit(counts)
    # fit_text = f"mu={mu:.1e} od={over_disp:.1e}"
    # if rsq:
    #     fit_text+=f" R2={rsq:.2e}"
    (obs, occ) = np.unique(counts, return_counts=True)
    x_max = obs[occ > 5].max()
    
    #mu_pois = poisson_fit(counts)

    x_fit = np.linspace(0, x_max, 500)
    y_fit_NB = nbinom.pmf(x_fit, n, p)
    #y_fit_pois = poisson.pmf(x_fit, mu_pois)
    
    traces = []
    
    traces.append(go.Histogram(
        histnorm =  'probability',
        x=counts,
        showlegend=False,
        xbins=dict(size=1)
    ))
    traces.append(go.Scatter(
        name= 'NB fitted',
        showlegend=False,
        x=x_fit,
        y=y_fit_NB
    ))
    # traces.append(go.Scatter(
    #     name= 'Poisson fitted',
    #     showlegend=False,
    #     x=x_fit,
    #     y=y_fit_pois,
    #     line=dict(dash='dashdot')
        
    # ))
    
    return traces
    
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir-count-matrices", type=Path, required=True)
    parser.add_argument("-o", type=Path, required=True)
    args = parser.parse_args()


    counts_dir = args.dir_count_matrices    #Path("data/simulated_data_large/kmer_quantification/count_matrices/")
    glob_str = "*.tsv"
    cluster_counts = counts_dir.glob(glob_str)
    count_files = list(x for x in cluster_counts if not x.name.startswith("counts_all"))
    outdir = args.o #Path("data/simulated_data_large/results/13_count_distribution")
    outdir.mkdir(parents=True, exist_ok=True)


    if count_files == []:
        raise FileExistsError(counts_dir/glob_str)

    for count_file in count_files:
        bgc_catalogue = count_file.stem
        print(bgc_catalogue, file=sys.stderr)

        df_count = pd.read_csv(count_file, sep="\t", index_col=0)
        df_plot = df_count.reset_index().rename(columns={'index':'kmer'}).melt(id_vars="kmer", var_name="sample", value_name="counts")
        df_plot[["dataset","sample"]] = df_plot["sample"].str.split(".", expand=True)

        for name_ds, df_ds in df_plot.groupby("dataset"):
            # perform within sample outlier detection.
            df_ds_outlier = pd.concat(
                df_sample.assign(is_outlier = is_outlier_iqr(df_sample.counts))
                for name, df_sample in df_ds.groupby("sample")
            )
            colors = ['#636efa', '#EF553B']

            fig_comp = make_subplots(
                rows=2, cols=2,
                specs=[
                    [{"colspan": 2}, None],
                    [{}, {}],
                    ],
                row_width = [0.4, 0.6],
            vertical_spacing=0.1,
                subplot_titles = ["Outlier overview", "Fitted NB with outliers", "Fitted NB without outliers"]
                #shared_yaxes=True,
            )
            
            full_box = go.Box(
                y = df_ds_outlier["sample"].values,
                x = df_ds_outlier["counts"].values,
                x0 = ' ',
                y0 = ' ',
                orientation="h",
                showlegend=False,
                line = {'color': 'rgba(150,150,150,1)'},
                marker = {'color': 'rgba(150,150,150,0)'}
            )
            fig_comp.add_trace(
                full_box,
                row=1, col=1
            )

            for is_outlier, df_plot_i in df_ds_outlier.groupby("is_outlier"):
                distrubution_box_plot = go.Box(
                    y = df_plot_i["sample"].values,
                    x = df_plot_i["counts"].values,
                    x0 = ' ',
                    y0 = ' ',
                    fillcolor = 'rgba(255,255,255,0)', #hide box
                    legendgroup = str(is_outlier),
                    name = str(is_outlier),
                    line = {'color': 'rgba(255,255,255,0)'}, #hide box
                    marker = {'color': colors[is_outlier],  'opacity':0.5},
                    offsetgroup = "points",
                    orientation = 'h',
                    pointpos = 0.5,
                    jitter=0.5,
                    alignmentgroup = 'True',
                    boxpoints = 'all',
                    #showlegend = cat_i == 0,
                    #boxmean='sd'
                )
                fig_comp.add_trace(
                    distrubution_box_plot,
                    row=1, col=1
                )

            q1, q3 = df_ds_outlier.counts.quantile([0.25, 0.75]).values.flatten()
            IQR = q3-q1
            max_x = max(5, q3+3*IQR)

            # Fit the NB on the full dataset
            for trace in generate_NB_trace(df_ds_outlier.counts):
                fig_comp.add_trace(
                    trace,
                    row=2, col=1
                )

            # Fit the NB on the non-outlier dataset
            for trace in generate_NB_trace(df_ds_outlier.query("not is_outlier").counts):
                fig_comp.add_trace(
                    trace,
                    row=2, col=2
                )   


            fig_comp.update_layout(xaxis2 = dict(range = (0, max_x)), xaxis3 = dict(range= (0, max_x)))
            fig_comp.update_layout(
                title = f"Count Distributions for BGC-group: {bgc_catalogue}<br><sup>{name_ds}",
                boxmode="group",
                legend = dict(title = "Is Outlier")
            )
            out_subdir = outdir/bgc_catalogue
            out_subdir.mkdir(parents=True, exist_ok=True)
            fig_comp.write_image(out_subdir/(name_ds+"_"+bgc_catalogue+".png"), height=600)#, height=800, width=600)
            scope._shutdown_kaleido() #attempt to fix kaleido hanging


