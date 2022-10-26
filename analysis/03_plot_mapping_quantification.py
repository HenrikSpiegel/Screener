#!/usr/bin/env python3
import glob
import pandas as pd
import os, sys
from typing import List
import numpy as np
import plotly.express as px

from qsub_modules.qsub_base import Base

def gather_reference_quantification(
        quantification_dir: str="data/simulated_data/quantification_map/", 
        expected_labels: List[str]=None):
    
    summary_files = glob.glob(os.path.join(quantification_dir, "*/cmseq_summation.tsv"))
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

if __name__ == "__main__":

    import configparser
    config = configparser.ConfigParser()
    config.read("config/project_config.ini")

    plot_width = 800
    plot_height = 500
    import pathlib
    input_dir = "data/simulated_data/quantification_map/"
    outdir = "results/03_mapping_quantification"
    fp_df       =   os.path.join(outdir, "map_quantification_summary.csv")
    fp_plot_avg =   os.path.join(outdir, "average_quantificatoon_errors.png")
    fp_plot_per =   os.path.join(outdir, "perBGC_quantification_errors.png")
    fp_plot_res =   os.path.join(outdir, "perBGC_residuals.png")

    pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)

    print("Running: " + __file__, "outdir -> "+outdir, sep="\n",file=sys.stderr)

    #The expected range of simulations: # TODO: add as project config.
    rstep = config.getfloat("Simulation", "ReadsGBStep")
    rmin = config.getfloat("Simulation", "ReadsGBMin")
    rmax = config.getfloat("Simulation", "ReadsGBMax")
    precision = config.get("Simulation", "ReadsGBStep").split(".")[1].__len__()
    gbs_to_run = np.round(np.arange(rmin, rmax+rstep, rstep), precision).tolist()

    if config.get("Simulation", "ReadsGBExtra", fallback=None):
        gbs_extra = [float(x) for x in config.get("Simulation", "ReadsGBExtra").split(",")]
        gbs_to_run.extend(gbs_extra)
    expected_labels = [Base.gen_prefix(gb) for gb in gbs_to_run]

    print("Collecting quantification files and calculating error metrics", file=sys.stderr)
    df_ref_quant = gather_reference_quantification(input_dir, expected_labels)
    df_summary = calculate_error_metrics(df_ref_quant)
    print("writting -> "+fp_df, file=sys.stderr)
    df_summary.to_csv(fp_df)

    print("Generating and writing plots", file=sys.stderr)
    df_avg = df_summary\
        .groupby("readsGB")[["AE","RAE"]].mean().reset_index()\
        .melt(id_vars="readsGB",value_vars=["AE","RAE"], var_name="Metric")
    fig_avg = px.line(df_avg, x="readsGB",y="value",line_dash="Metric",
              labels={"value": "Error", "readsGB": "Sample size (Gb)"},
              markers=True,
             title="Average quantification errors<br><sup>Quantification By mapping<sup>")
    fig_avg.update_layout(xaxis_range=[0,0.35], yaxis_range=[0,4])
    print("writting -> "+fp_plot_avg, file=sys.stderr)
    fig_avg.write_image(fp_plot_avg, width=plot_width, height=plot_height)

    fig_pBGC = px.line(df_summary.sort_values("readsGB"), x="readsGB",y="RAE", color="Contig",
        labels={"value": "Error", "readsGB": "Sample size (Gb)","Contig":"BGC"},
              markers=True,
             title="Per BGC quantification errors<br><sup>Quantification By mapping<sup>")
    fig_pBGC.update_layout(xaxis_range=[0,0.35], yaxis_range=[0,4])
    print("writting -> "+fp_plot_per, file=sys.stderr)
    fig_pBGC.write_image(fp_plot_per, width=plot_width, height=plot_height)

    fig_res = px.scatter(df_summary, x="readsGB",y="std_res",
            labels={"std_res": "Standardized residuals", "readsGB": "Sample size (Gb)"},
            title="Residual Plot for each BGC<br><sup>Difference between Median Depth and expected BGC count<sup>"
          )
    print("writting -> "+fp_plot_res, file=sys.stderr)
    fig_res.write_image(fp_plot_res, width=plot_width, height=plot_height)