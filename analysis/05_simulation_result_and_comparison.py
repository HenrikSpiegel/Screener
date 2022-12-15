
import pandas as pd

import plotly.express as px
from pathlib import Path
import json
from plotly.subplots import make_subplots
import plotly.graph_objects as go

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--simulation-overview", required=True, type=Path)
    parser.add_argument("--total-abundances", required=True, type=Path)
    parser.add_argument("--family-dump", required=True, type=Path)
    parser.add_argument("-o", required=True, type=Path)
    args = parser.parse_args()

    fp_simulation = args.simulation_overview  #Path("data/simulated_data_large/camisim/simulation_overview_full.tsv")
    fp_abundances = args.total_abundances     #Path("data/simulated_data_large/abundances/total_abundances.csv")
    fp_json_families = args.family_dump       #Path("data/simulated_data_large/mcl_clustering/out.blast_result.mci.I40.json")
    outdir = args.o                           #Path("data/simulated_data_large/results/05_simulation_result_and_comparison")
    outdir.mkdir(parents=True, exist_ok=True)

    ### Combine simulation overview and catalogue groupings to get expected values.
    df_simulation = pd.read_csv(fp_simulation, sep="\t")

    dict_families = json.loads(fp_json_families.read_text())
    by_sample_grouped = df_simulation.groupby(["dataset", "sample"])
    rows = []
    for name, df_g in by_sample_grouped:
        group_rows = [
            name + (cat, df_g.loc[df_g.ncbi.isin([member.rsplit(".",1)[0] for member in cat_members]),'expected_average_coverage'].sum())
            for cat, cat_members in dict_families.items()]
        rows.extend(group_rows)
    df_catalogue_expect = pd.DataFrame(rows, columns = ["dataset", "sample","catalogue_name", "expected_average_coverage"])

    ### Load abundances.
    df_abundances = pd.read_csv(fp_abundances)
    df_abundances[["dataset", "sample"]] = df_abundances.dataset_sample.str.split(".", expand=True)

    ### Combine and calculate errors:
    df_combined = df_catalogue_expect.merge(df_abundances, on = ["dataset", "sample","catalogue_name"])
    df_err_long = df_combined\
        .melt(
            id_vars = ["catalogue_name","method", "dataset","sample", "expected_average_coverage"],
            value_vars = ["negbinom_mu", "median", "negbinom_mu_corr", "median_corr"], 
            value_name = "estimate",
            var_name = "estimate_agg")
    df_err_long["RE"] = (df_err_long.estimate-df_err_long.expected_average_coverage) / df_err_long.expected_average_coverage
    df_err_long["RAE"] = df_err_long["RE"].abs()

    #########################################################
    ################## Generate RE plot:#####################
    px.colors.qualitative.Plotly
    color_map = {
        'init':px.colors.qualitative.Plotly[0],
        'best':px.colors.qualitative.Plotly[1],
        'raw_5000':px.colors.qualitative.Plotly[2]
    }
    df_plot = df_err_long.query("estimate_agg == 'negbinom_mu_corr'")

    catagories, df_catagories = zip(*list(df_plot.groupby("catalogue_name")))
    titles = []
    for cat in catagories:
        titles.extend([cat, "Summarised"])
    titles.extend(["Combined catagories", "Total"])
    # Create figure
    fig_comp = make_subplots(
        rows=len(catagories)+1, cols=2,
        shared_yaxes=True,
        column_widths = [0.7, 0.3],
        subplot_titles=titles,
        #vertical_spacing=0.06,
        #horizontal_spacing=0.03
    )


    for cat_i, (cat_name, df_group_cat) in enumerate(zip(catagories, df_catagories)):
        
        for method, df_method in df_group_cat.groupby("method"):
            df_method.sort_values("dataset", inplace=True)
            
            fig_comp.add_trace(
                go.Box(
                    x = df_method["dataset"].values,
                    y = df_method["RE"].values,
                    fillcolor = 'rgba(255,255,255,0)', #hide box
                    legendgroup = method,
                    name = method,
                    line = {'color': 'rgba(255,255,255,0)'}, #hide box
                    marker = {'color': color_map[method]},
                    offsetgroup = method,
                    orientation = 'v',
                    pointpos = 0,
                    jitter=0.3,
                    alignmentgroup = 'True',
                    boxpoints = 'all',
                    showlegend = cat_i == 0,
                    boxmean='sd'
                ),
                row=cat_i+1, col=1
            )
            
            fig_comp.add_trace(
                go.Box(
                    x = ["Across datasets" for i in range(len(df_method))],
                    y = df_method["RE"].values,
                    fillcolor = 'rgba(255,255,255,0)', #hide box
                    legendgroup = method,
                    name = method,
                    line = {'color': 'rgba(255,255,255,0)'}, #hide box
                    marker = {'color': color_map[method]},
                    offsetgroup = method,
                    orientation = 'v',
                    pointpos = 0,
                    jitter=0.5,
                    alignmentgroup = True,
                    boxpoints = 'all',
                    showlegend = False,
                    #boxmean='sd'
                ),
                row=cat_i+1, col=2
            )
    # add final summarising row.
    for method, df_method in df_plot.groupby("method"):
            df_method.sort_values("dataset", inplace=True)
            
            fig_comp.add_trace(
                go.Box(
                    x = df_method["dataset"].values,
                    y = df_method["RE"].values,
                    fillcolor = 'rgba(255,255,255,0)', #hide box
                    legendgroup = method,
                    name = method,
                    line = {'color': 'rgba(255,255,255,0)'}, #hide box
                    marker = {'color': color_map[method]},
                    offsetgroup = method,
                    orientation = 'v',
                    pointpos = 0,
                    jitter=0.3,
                    alignmentgroup = 'True',
                    boxpoints = 'all',
                    showlegend = cat_i == 0,
                    boxmean='sd'
                ),
                row=cat_i+2, col=1
            )
            
            fig_comp.add_trace(
                go.Box(
                    x = ["Across datasets" for i in range(len(df_method))],
                    y = df_method["RE"].values,
                    fillcolor = 'rgba(255,255,255,0)', #hide box
                    legendgroup = method,
                    name = method,
                    line = {'color': 'rgba(255,255,255,0)'}, #hide box
                    marker = {'color': color_map[method]},
                    offsetgroup = method,
                    orientation = 'v',
                    pointpos = 0,
                    jitter=0.5,
                    alignmentgroup = True,
                    boxpoints = 'all',
                    showlegend = False,
                    #boxmean='sd'
                ),
                row=cat_i+2, col=2
            )
    fp_fig_comp = outdir/"RE_overview.png"
    fig_comp.update_layout(
        height=150*len(catagories)+1, 
        boxmode='group', 
        title="Overview of Relative error for each cluster across datasets.")
    fig_comp.write_image(fp_fig_comp,scale=2)

    ###########################################
    ########### DISTRIBUTION PLOTS ############
    fig_ecov_genomes = px.strip(
        df_simulation,
        x="dataset",
        y="expected_average_coverage",
        color="ncbi",
        log_y=True, 
        labels = {
            'expected_average_coverage':'Expected Average Coverage',
            'dataset':'Dataset'
        }
    )
    fig_ecov_genomes.update_layout(
        title="Distribution of expected coverage of the genomes over the datasets",
        showlegend=False)
    fig_ecov_genomes.write_image(outdir/"ecov_genomes.png", height=550, width=1000)
    
    ###

    fig_dist_genomes = px.strip(
        df_simulation,
        x="dataset",
        y="distribution",
        color="ncbi", 
        labels={
            'dataset':'Dataset',
            'distribution':'Relative Frequency'
        }
    )
    fig_dist_genomes.update_layout(
        title="Distribution of the relative frequencies of the genomes over the datasets",
        showlegend=False)
    fig_dist_genomes.write_image(outdir/"dist_genomes.png", height=550, width=1000)

    ###
    fig_ecov_cat = px.strip(
        df_catalogue_expect,
        x="dataset",
        y="expected_average_coverage",
        color="catalogue_name",
        log_y=True,
        labels={'expected_average_coverage': 'Expected Average Coverage',
                'catalogue_name':'Catalogue Name',
                'dataset':'Dataset'
            }
    )
    fig_ecov_cat.update_layout(
        title="Distribution of expected coverages for the various catalogues over the datasets.",
        showlegend=True,
        width=900,
        height=650
    )
    fig_ecov_cat.write_image(outdir/"ecov_catalogues.png", height=650, width=1000)

