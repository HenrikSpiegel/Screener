
import pandas as pd
import sys
import numpy as np
import plotly.express as px
from pathlib import Path


#NOTE: Teribly hardcoded paths :/

import configparser
def get_simulation_overview():
    dfs = []
    for summary_file in Path("data/simulated_data/camisim/").glob("*GB/simulation_overview.csv"):
        df = pd.read_csv(summary_file, sep=",")
        df["dataset"] = summary_file.parent.name
        dfs.append(df)
    df_simulation = pd.concat(dfs)
    
    #Grab config data:
    camisim_config = configparser.ConfigParser()
    #Here we assume only size (ie. readdepth) changes.
    config_file = Path("data/simulated_data/camisim/configs/").glob("*_config.ini").__next__()
    camisim_config.read(config_file)
    df_simulation["err_type"]=camisim_config.get("ReadSimulator","type")
    try:
        df_simulation["err_profile"]=camisim_config.getfloat("ReadSimulator","profile")
    except:
        df_simulation["err_profile"]=camisim_config.get("ReadSimulator","profile")
    df_simulation["fragment_size"]=camisim_config.getint("ReadSimulator","fragments_size_mean")
    
    return df_simulation

# get quantifications.

def get_kmer_summation():
    dfs = []
    for summary_file in Path("data/simulated_data/quantification_kmer/").glob("*GB/sample_*/kmer_summation.tsv"):
        try:
            dataset = summary_file.parent.parent.name
            sample  = summary_file.parent.name
            df = pd.read_csv(summary_file, sep="\t")
            df["sample"] = sample
            df["dataset"] = dataset
            df["quantification"] = "kmerCount"
            dfs.append(df)
        except Exception as err:
            print(summary_file)
            raise err
    df_summation = pd.concat(dfs)
    return df_summation

def get_map_summation():
    dfs = []
    for summary_file in Path("data/simulated_data/quantification_map/").glob("*GB/sample_*/cmseq_summation.tsv"):
        try:
            dataset = summary_file.parent.parent.name
            sample  = summary_file.parent.name
            df = pd.read_csv(summary_file, sep="\t")
            df["sample"] = sample
            df["dataset"] = dataset
            df["quantification"] = "mapping(filtered)"
            dfs.append(df)
        except Exception as err:
            print(summary_file)
            raise err
    for summary_file in Path("data/simulated_data/quantification_map/").glob("*GB/sample_*/cmseq_summation_raw.tsv"):
        try:
            dataset = summary_file.parent.parent.name
            sample  = summary_file.parent.name
            df = pd.read_csv(summary_file, sep="\t")
            df["sample"] = sample
            df["dataset"] = dataset
            df["quantification"] = "mapping(unfiltered)"
            dfs.append(df)
        except Exception as err:
            print(summary_file)
            raise err
    df_summation = pd.concat(dfs)
    return df_summation

def get_name_for_run():
    camisim_config = configparser.ConfigParser()
    #Here we assume only size (ie. readdepth) changes.
    config_file = Path("data/simulated_data/camisim/configs/").glob("*_config.ini").__next__()
    camisim_config.read(config_file)
    
    sim       = camisim_config.get("ReadSimulator","type")
    profile   = camisim_config.get("ReadSimulator","profile").replace(".","")
    frag_size = camisim_config.get("ReadSimulator","fragments_size_mean")
    #frag_size = frag_size
    seed      = camisim_config.get("Main","seed")
    sample_type = camisim_config.get("community0","mode")
    name = f"{sim}_{profile}_frag_{frag_size}_seed_{seed}_sample_{sample_type}"
    return name

if __name__ == '__main__':
    print("Running: " + __file__, file=sys.stderr)

    # Collect data
    df_sum = pd.concat([get_kmer_summation(), get_map_summation()])
    df_sum[["ncbi","region"]] = df_sum.Contig.str.rsplit(".",1, expand=True)
    df_sum.ncbi = df_sum.ncbi.str.replace("\.","_")

    df_sim = get_simulation_overview()
    df_sim.rename(columns={"samle":"sample"}, inplace=True)

    df_comb = pd.merge(df_sum, df_sim, how="inner", on=["dataset","sample","ncbi"])
    df_comb["readsGB"] = df_comb.dataset.str.replace("GB","")
    df_comb["readsGB"] = df_comb.readsGB.str.replace("_",".").astype(float)
    df_comb.sort_values(["readsGB","Contig"], inplace=True)

    runname = get_name_for_run()
    res_dir = Path("results/") / Path(__file__).stem / runname
    res_dir.mkdir(parents=True, exist_ok=True)
    print(f"outdir {res_dir}", file=sys.stderr)

    df_comb.to_csv(res_dir / "combined_simulation_and_sumation.csv")

    # Simulation overview plots
    fig = px.box(df_comb, x="readsGB", y="distribution", color="ncbi", 
       title="Distribution of genomes<br><sup>Over the samples<sup>",
        labels={'ncbi':'Genome', 'readsGB':'Sample size (Gb)'},
      points="all")
    fig.write_image(str(res_dir/"genome_distribution.png"), height=500, width=800, scale=3)
    fig

    fig = px.box(df_comb, x="readsGB", y="expected_average_coverage", color="ncbi",
       title="Average coverage of genomes<br><sup>Over the samples<sup>",
          labels={'ncbi':'Genome', 'readsGB':'Sample size (Gb)',
                  'expected_average_coverage':'Average Coverage'
                 },
      points="all")
    fig.write_image(str(res_dir/"genome_coverage.png"), height=500, width=800, scale=3)


    #Wrangle the data so we can get the raw count and corrected count from the kmerquant rows.
    df_plot = df_comb.melt(
        id_vars=["Contig","expected_average_coverage","readsGB", "quantification"], 
        value_vars=["Depth median", 'Depth Error Corrected', 'Depth Error/Edge Corrected'],
        value_name="estimate",
        var_name="sub_quantification"
    )
    df_plot.dropna(subset=["estimate"], inplace=True)
    df_plot.reset_index(drop=True, inplace=True)
    quantification_name = [q if q.startswith("mapping") else f"{q}({v})" for q,v in zip(df_plot["quantification"], df_plot["sub_quantification"])]
    df_plot["quantification_name"] = quantification_name
    df_plot.sort_values(["Contig", "quantification_name"], inplace=True)

    # add error values
    df_plot["error"] = df_plot["estimate"] - df_plot["expected_average_coverage"]
    df_plot["RE"]    = df_plot["error"]/df_plot["expected_average_coverage"]



    # Error comparison plots:
    fig = px.box(df_plot, x="Contig", color='quantification_name', y="error",facet_row="readsGB", 
        height=2000,width=1000, points="all",
        title="Distribution of errors for each contig.<br><sup>For multiple samples<sup>",
        labels={'error':'Estimate - Expected', 'quantification_name':'Method'}
        )
    fig.update_yaxes(matches=None)
    fig.write_image(str(res_dir/"error_plot.png"), height=2000, width=1000, scale=5)

    ffig = px.box(df_plot, x="Contig", color='quantification_name', y="RE",facet_row="readsGB", 
       height=2000,width=1000, points="all",
       title="Distribution of relative error for each contig.<br><sup>For multiple samples<sup>",
       labels={'RE':'(Estimate - Expected)/Expected', 'quantification_name':'Method'}
      )
    fig.write_image(str(res_dir/"relative_error_plot.png"), height=2000, width=1000, scale=5)

    print("-- finished "+__file__, file=sys.stderr)
    