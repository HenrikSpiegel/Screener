import glob
import pathlib
import pandas as pd
import os, sys
from typing import List, Literal
import numpy as np
import plotly.express as px
import plotly.graph_objects as go




def gather_quantification(
        fuzzy_filename: Literal[
            "data/simulated_data/quantification_map/*/cmseq_summation.tsv",
            "data/simulated_data/quantification_map/*/kmer_summation.tsv"],
        expected_labels: List[str]=None,
):
    
    summary_files = glob.glob(fuzzy_filename)
    summary_labels = [x.rsplit("/",2)[1] for x in summary_files]
    
    if expected_labels:
        if not set(expected_labels) == set(summary_labels):
            err_msg = f"\
expected summary files != found summary files.\n\
Expected labels({len(expected_labels)}):{expected_labels}\n\
Found labels({len(summary_labels)}):{summary_labels}\n\
Missing expected: {set(expected_labels).difference(set(summary_labels))}"
            raise ValueError(err_msg)
    
    dataframes = []
    for file, label in zip(summary_files, summary_labels):
        readsGB = float(label[:-2].replace("_", "."))
        df = pd.read_csv(file, sep="\t")
        df["readsGB"] = readsGB
        dataframes.append(df)        
    print(f"Gathered summary from ({len(summary_labels)}): {summary_labels}", file=sys.stderr)
    return pd.concat(dataframes)

def expected_count_eq_ratio(readsGB:List[int], input_genomes_fasta: str = "data/simulated_data/input_genomes/combined.fa"):
    total_length = 0
    for line in open(input_genomes_fasta, "rb"):
        if line[0] == 62: #bit for '>'
            #print(total_length)
            #print(line)
            continue
        total_length += len(line.strip())
    print(f"Combined genome size ({total_length})bp", file=sys.stderr)
    expected_counts = [(x*10**9)/total_length for x in readsGB]
    return expected_counts


def calculate_error_metrics(df: pd.DataFrame) -> pd.DataFrame:
    df_out = df.copy()
    
    fp_inputcombined = "data/simulated_data/input_genomes/combined.fa"
    df_out["expected_value"] = expected_count_eq_ratio(df["readsGB"], fp_inputcombined)
    df_out["expected_value_int"] = df_out["expected_value"].round(0)
    
    df_out["err"] = df_out["Depth median"] - df_out["expected_value"]
    df_out["AE"] = np.abs(df_out["err"])
    df_out["RAE"] = df_out["AE"]/df_out["expected_value"] 
    
    #standardized residuals
    df_out["std_res"] = df_out["err"] / np.std(df_out["err"], ddof=2)
    return df_out

def gen_prefix(reads_gb:int) -> str:
    return str(reads_gb).replace(".","_")+"GB"


if __name__ == "__main__":
    sys.stderr.write(f"---Running: {__file__}\n")
    import configparser
    config = configparser.ConfigParser()
    config.read("config/project_config.ini")

    sys.stderr.write("Prefacing...\n")
    rstep = config.getfloat("Simulation", "ReadsGBStep")
    rmin = config.getfloat("Simulation", "ReadsGBMin")
    rmax = config.getfloat("Simulation", "ReadsGBMax")
    precision = config.get("Simulation", "ReadsGBStep").split(".")[1].__len__()
    gbs_to_run = np.round(np.arange(rmin, rmax+rstep, rstep), precision).tolist()

    if config.get("Simulation", "ReadsGBExtra", fallback=None):
        gbs_extra = [float(x) for x in config.get("Simulation", "ReadsGBExtra").split(",")]
        gbs_to_run.extend(gbs_extra)
        
    expected_labels = [gen_prefix(x) for x in gbs_to_run]

    output_dir = "results/04_comparison_of_initial_methods"
    pathlib.Path(output_dir).mkdir(parents=True, exist_ok=True)
    fp_csv_raw      = os.path.join(output_dir, "raw_stats.csv")
    fp_csv_summary  = os.path.join(output_dir, "summarised_stats.csv")
    fp_fig_count    = os.path.join(output_dir, "initial_compare_count.png")
    fp_fig_ae       = os.path.join(output_dir, "initial_compare_ae.png")
    fp_fig_rae      = os.path.join(output_dir, "initial_compare_rea.png")
    plot_width = 800
    plot_height = 500

    sys.stderr.write("Gathering data...\n")
    df_map = gather_quantification("data/simulated_data/quantification_map/*/cmseq_summation.tsv", expected_labels)
    df_map_raw = gather_quantification("data/simulated_data/quantification_map/*/cmseq_summation_raw.tsv", expected_labels)
    df_kmer = gather_quantification("data/simulated_data/quantification_kmer/*/kmer_summation.tsv",expected_labels)

    df_map["type"] = "mapping_MapQ30"
    df_map_raw["type"] = "mapping_MapQ00"
    df_kmer["type"] = "K-mer"

    df_comb = pd.concat([df_map, df_kmer, df_map_raw])
    df_comb_stat = calculate_error_metrics(df_comb)

    df_comb_summary = df_comb_stat.groupby(["type","readsGB"])[["expected_value", "Depth median","AE", "RAE"]].mean().reset_index()
    df_comb_summary_long = df_comb_summary.melt(id_vars=["type","readsGB"], value_vars=["expected_value","Depth median","AE","RAE"])

    df_comb_stat.to_csv(fp_csv_raw, index=False)
    df_comb_summary.to_csv(fp_csv_summary, index=False)

    sys.stderr.write(f"Generating and writing plots -> {output_dir} ...\n")

    df_expected = df_comb_summary_long.query("variable == 'expected_value'")
    trace_expected = go.Scatter(
        x=df_expected["readsGB"],
        y=df_expected["value"],
        mode='lines',
        name='Expected',
        line=dict(color='black', width=2,
                                dash='5px'))
    
    fig_count = px.line(df_comb_summary_long.query("variable == 'Depth median'"), 
        x="readsGB", y="value", color="type",
        labels={'value':'Cluster Count', 'readsGB':'Sample size (Gb)', 'type':'Method'},
        title="Comparison of quantification by initial methods.<br><sup>Cluster counts<sup>")
    fig_count.add_trace(trace_expected)
    fig_count.write_image(fp_fig_count, height=plot_height, width=plot_width)

    fig_rae = px.line(df_comb_summary_long.query("variable == 'RAE'"), 
        x="readsGB", y="value", color="type",
       labels={'value':'Relative Absolute Error (RAE)', 'readsGB':'Sample size (Gb)', 'type':'Method'},
       title="Comparison of quantification by initial methods.<br><sup>Relative absolute error of quantification<sup>")
    fig_rae.write_image(fp_fig_rae, height=plot_height, width=plot_width)

    fig_ae = px.line(df_comb_summary_long.query("variable == 'AE'"), 
        x="readsGB", y="value", color="type",
       labels={'value':'Absolute Error (AE)', 'readsGB':'Sample size (Gb)', 'type':'Method'},
       title="Comparison of quantification by initial methods.<br><sup>Absolute error of quantification<sup>")
    fig_ae.write_image(fp_fig_ae, height=plot_height, width=plot_width)
    sys.stderr.write(f"Finished -> {__file__}\n")