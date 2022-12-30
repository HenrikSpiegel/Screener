from pathlib import Path
import pandas as pd
import json
from Bio import SeqIO

import plotly.express as px
import plotly.graph_objects as go
import plotly
import numpy as np

def get_map_quant(mapping_top_dir: Path):
    map_quant_glob = mapping_top_dir.glob("*GB/sample*/cmseq_summation.tsv")
    summation_frames = []
    for sum_file in map_quant_glob:
        dataset_name = sum_file.parent.parent.name
        sample_name = sum_file.parent.name
        summation_frames.append(
            pd.read_csv(sum_file, sep="\t")\
            .assign(dataset=dataset_name, sample=sample_name)
        )
    df = pd.concat(summation_frames).reset_index(drop=True)

    df.rename(columns = {
        'Contig': 'catalogue_name',
        'Depth avg': 'depth_avg',
        'Depth median': 'depth_median'
    }, inplace=True)
    df_out = df.loc[:, ["catalogue_name", "dataset","sample", "depth_avg", "depth_median"]]
    df_out['method'] = "MapToRef"
    return df_out

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument("--dir-map-quant", required=True, type=Path)
    parser.add_argument("--kmer-abundances", required=True, type=Path)
    parser.add_argument("--camisim-overview", required=True, type=Path)
    parser.add_argument("--cluster-json", required=True, type=Path)
    parser.add_argument("-o", type=Path, required=True)
    args = parser.parse_args()

    map_quant_dir = args.dir_map_quant          #Path("../data/simulated_data_large/map_quantification/")
    kmer_quant_fp = args.kmer_abundances         #Path("../data/simulated_data_large/abundances/total_abundances.csv")
    camisim_overview = args.camisim_overview     #Path("../data/simulated_data_large/camisim/simulation_overview_full.tsv")
    cluster_def_json = args.cluster_json         #Path("../data/simulated_data_large/mcl_clustering/out.blast_result.mci.I40.json")
    outdir = args.o

    #Mapping quantification
    df_mapping = get_map_quant(map_quant_dir)

    #Kmer quantification:
    df_kmer = pd.read_csv(kmer_quant_fp)
    df_kmer.loc[:,"method"] = "kmer_"+df_kmer.loc[:,"method"]
    df_kmer[["dataset", "sample"]] = df_kmer.dataset_sample.str.split(".", expand=True)

    

    #Extract relevant columns and combine.
    keep_cols = ["catalogue_name","dataset","sample", "method", "estimate"]
    keep_methods = ["MapToRef","kmer_random","kmer_init","kmer_best"]
    
    df_mapping_stack = df_mapping.rename(columns = {'depth_median':'estimate'})[keep_cols]
    df_kmer_stack = df_kmer.rename(columns = {'negbinom_mu_corr':'estimate'})[keep_cols]
    df_combined = pd.concat([df_mapping_stack, df_kmer_stack])

    df_estimates = df_combined.query(f"method in {keep_methods}")

    #Get per catalogue entry expected abundances:
    df_simulation = pd.read_csv(camisim_overview, sep="\t")
    dict_families = json.loads(cluster_def_json.read_text())

    by_sample_grouped = df_simulation.groupby(["dataset", "sample"])
    rows = []
    for name, df_g in by_sample_grouped:
        group_rows = [
            name + (cat, df_g.loc[df_g.ncbi.isin([member.rsplit(".",1)[0] for member in cat_members]),'expected_average_coverage'].sum())
            for cat, cat_members in dict_families.items()]
        rows.extend(group_rows)
    df_catalogue_expect = pd.DataFrame(rows, columns = ["dataset", "sample","catalogue_name", "expected_average_coverage"])

    #Combine expected and estimates and calculate errors:
    df_total = df_estimates.merge(df_catalogue_expect, on = ["dataset","sample","catalogue_name"]).reset_index(drop=True)
    df_total["RE"] = (df_total.estimate-df_total.expected_average_coverage) / df_total.expected_average_coverage
    df_total["RAE"] = df_total["RE"].abs()
    df_total["dataset_size"] = [float(x.strip("GB").replace("_",".")) for x in df_total.dataset]

    color_map = { #Yep this is wierd setup
        method: ["blue", "red","goldenrod", "green"][i]
        for i, method in enumerate(["MapToRef","kmer_random","kmer_init","kmer_best"])
        }

    
    #### Generating plots

    # RE plots:
    outdir_re = outdir/"re_plots"
    outdir_re.mkdir(parents=True, exist_ok=True)


    for catalogue in df_total.catalogue_name.drop_duplicates():
        df_plot = df_total.query(f"catalogue_name == '{catalogue}'")
        df_plot_sum = df_plot.groupby(["dataset", "dataset_size","method"])["RE"].agg([np.mean, np.std]).reset_index()
        fig = px.line(
            df_plot_sum, 
            x="dataset_size",
            log_x=True,
            y="mean",
            error_y="std",
            color="method",
            color_discrete_map=color_map,
            title=f"Catalogue: {catalogue} - Comparison of quantification methods<br><sup>Mean relative error (MRE) and the standard deviation of the MRE",
            labels={
                'mean': 'Mean Relative Error',
                'dataset_size' :"Dataset Size in GB"
            }
        )
        fig.update_yaxes(range=[-1, 2])
        outfile = outdir_re / (catalogue+".png")
        fig.write_image(outfile, scale=1.5)

    #est / obs plots
    outdir_est = outdir/"estimate_plots"
    outdir_est.mkdir(parents=True, exist_ok=True)


    for catalogue in df_total.catalogue_name.drop_duplicates():
        df_plot = df_total.query(f"catalogue_name == '{catalogue}'")
        df_plot_sum = df_plot.groupby(["dataset","dataset_size","method"])["estimate"].agg([np.mean, np.std]).reset_index()
        fig = px.line(
            df_plot_sum, 
            x="dataset_size",
            log_x=True,
            log_y=True,
            y="mean",
            error_y="std",
            color="method",
            color_discrete_map=color_map,
            title=f"Catalogue: {catalogue} - Comparison of quantification methods<br><sup>Mean estimates across samples and standard deviation.",
            labels={
                'mean': 'Mean Estimate',
                'dataset_size': "Dataset Size in GB"
            }
        )
        fig.update_traces(opacity=0.5)
        
        df_expect = df_plot.groupby(["dataset","dataset_size","method"])["expected_average_coverage"].agg([np.mean, np.std]).reset_index()
        fig.add_trace(
            go.Scatter(
                x=df_expect["dataset_size"].values,
                y=df_expect["mean"].values,
                name="Expected",
                line = dict(color='black', width=2, dash='dash')
            )
        )
        #fig.update_yaxes(range=[-1, 2])
        outfile = outdir_est / (catalogue+".png")
        fig.write_image(outfile, scale=1.5)