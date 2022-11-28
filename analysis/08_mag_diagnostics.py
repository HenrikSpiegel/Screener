import sys
import numpy as np
import pandas as pd
from pathlib import Path
import configparser
import json
from tqdm import tqdm
from typing import Sequence
import plotly.figure_factory as ff

#from scripts.kmer_summarise import error_correction, edgeloss_correction

def error_correction(values: Sequence[float], k:int=21, error_rate:float=0.03):
    correction_value = (1-error_rate)**k
    corrected_value = [x/correction_value for x in values]
    return corrected_value


def edgeloss_correction(values: Sequence[float], k:int=21, L:int=150):
    correction_value = 1 - ((k-1)/L)
    corrected_value = [x/correction_value for x in values]
    return corrected_value


def get_simulation_overview(file_simulation_overview_fuzzy="data/simulated_data/camisim/*GB/simulation_overview.csv"):
    dfs = []
    for summary_file in Path().glob(file_simulation_overview_fuzzy):
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

if __name__ == '__main__':
    
    # Getting expected values:
    df_simulation = get_simulation_overview("data/simulated_data/camisim/*GB/simulation_overview.csv")
    df_simulation["readsGB"]  = df_simulation.dataset.str.replace("_",".").str.rstrip("GB").astype(float)
    
    file_catalogue_grouping = Path("data/simulated_data/catalogues/family_dump.json")
    catalogue_groupings     = json.loads(file_catalogue_grouping.read_text())

    by_sample_grouped = df_simulation.groupby(["dataset", "sample"])
    rows = []
    for name, df_g in by_sample_grouped:
        group_rows = [
            name + (cat, df_g.loc[df_g.ncbi.isin([member.rsplit(".",1)[0].replace(".","_") for member in cat_members]),'expected_average_coverage'].sum())
            for cat, cat_members in catalogue_groupings.items()]
        rows.extend(group_rows)
    df_catalogue_expect = pd.DataFrame(rows, columns = ["dataset", "sample","catalogue_name", "expected_average_coverage"])
    
    # Collecting counts:
    dir_count_matrices = Path("data/simulated_data/kmer_quantification/count_matrices/")
    catalouge_count_files = [file for file in dir_count_matrices.glob("*.tsv") if not file.stem=="counts_all"]

    df_count = pd.concat(
        pd.read_csv(file, sep="\t", index_col=0)\
            .reset_index()\
            .rename(columns={'index':'kmer'})\
            .melt(id_vars=["kmer"], var_name='dataset_sample', value_name='count')\
            .assign(catalogue_name = file.stem)
        for file in catalouge_count_files
    )
    df_count[["dataset","sample"]] = df_count['dataset_sample'].str.split(".", expand=True)

    # count correction

    def combined_correction(values, k, L, error_rate):
        err_corrected = error_correction(values=values, k=k, error_rate=error_rate)
        
        return edgeloss_correction(
            values=err_corrected,
            k=k, L=L
        )
    k = 21 #load from config ...
    L = 148.5
    error_rate = 0.015
    df_count["count_corrected"] = combined_correction(df_count["count"], k, L, error_rate)


    # Merge expected and conut and add error values
    df_count_combined = pd.merge(df_count, df_catalogue_expect, on = ['dataset', 'sample', 'catalogue_name'])
    df_count_combined["RE"] = (df_count_combined['count_corrected'] - df_count_combined['expected_average_coverage']) / df_count_combined['expected_average_coverage']

   

    # collect MAGinator gene sets
    files_mag_gene_sets = Path("data/simulated_data/MAGinator/screened_flat/").glob("*_kmers.csv")
    df_genes_sets = pd.concat(
        pd.read_csv(file)\
            .assign(catalogue_name = file.stem.rsplit("_kmers",1)[0])
        for file in files_mag_gene_sets
    )

    # Generate visualizations
    catalogues = sorted(df_genes_sets.loc[:,'catalogue_name'].drop_duplicates())
    outdir = Path("results") / (Path(__file__).stem+"_05GB")
    #TODO: TEMP!!!
    print("WARNING: TEMPORARY extractions", file=sys.stderr)
    df_count_combined = df_count_combined[df_count_combined["dataset"] == "0_5GB"]
    outdir.mkdir(parents=True, exist_ok=True)

    for catalogue in (pbar := tqdm(catalogues)):
        pbar.set_description(f"[{catalogue}]: Processing data.")
        df_gene_i = df_genes_sets.loc[df_genes_sets['catalogue_name']==catalogue,:]
        df_i = df_count_combined.loc[df_count_combined['catalogue_name']==catalogue, ['kmer', 'dataset_sample', 'RE']]
        
        hist_data = [
            df_i.loc[:,'RE'].values,
            df_i.loc[df_i['kmer'].isin(df_gene_i['init']),'RE'].values,
            df_i.loc[df_i['kmer'].isin(df_gene_i['best']),'RE'].values
        ]

        MRAE = [np.mean(np.abs(x)) for x in hist_data]

        group_labels = [
            f'all_{len(df_i.kmer.drop_duplicates())}',
            'MAG_init',
            'MAG_best'
        ]

        mse_string = "MRAE: "+ " ".join([f"{label}: {value:0.1f}" for label, value in zip(group_labels, MRAE) ])

        pbar.set_description(f"[{catalogue}]: Processing plot.")
        # Create distplot with custom bin_size
        fig = ff.create_distplot(hist_data, group_labels, bin_size=.1,)
        fig.update_layout(
            title=f"Histogram with KDE curve showing distribution of relative errors for catelogue groupings.<br><sup>{catalogue} - {mse_string}",
            yaxis_title="Probability Density",
            xaxis_title="Relative Error",
            yaxis1=dict(range=[0,2])
            )
        fig.update_traces(opacity=0.75)

        if catalogue == "NZ_CP053893.1.region002":
            fig.update_xaxes(range=[-1,25])
        else:
            fig.update_xaxes(range=[-1,5])

        file_image = outdir / (catalogue+".png")
        fig.write_image(file_image, width=1000, height=800)
