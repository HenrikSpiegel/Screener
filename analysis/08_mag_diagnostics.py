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


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--simulation-overview", type=Path, required=True)
    parser.add_argument("--family-json", type=Path, required=True)
    parser.add_argument("--count-mats", type=Path, required=True)
    parser.add_argument("--mag-flat", type=Path, required=True)
    parser.add_argument("-o", type=Path, required=True)
    args = parser.parse_args()


    fp_simulation_overview = args.simulation_overview   #"data/simulated_data_init/camisim/*GB/simulation_overview.csv"
    bgc_family_dump = args.family_json          #Path("data/simulated_data_init/catalogues/family_dump.json")
    dir_count_matrices = args.count_mats        #Path("data/simulated_data_init/kmer_quantification/count_matrices/")
    dir_mag_flat = args.mag_flat        #Path("data/simulated_data_init/MAGinator/screened_flat/")

    outdir = args.o #("data/simulated_data_init/results/08_mag_diagnostics")

    # Getting expected values:
    df_simulation = pd.read_csv(fp_simulation_overview, sep="\t")
    df_simulation["readsGB"]  = df_simulation.dataset.str.replace("_",".").str.rstrip("GB").astype(float)
    
    file_catalogue_grouping = Path(bgc_family_dump)
    catalogue_groupings     = json.loads(file_catalogue_grouping.read_text())

    by_sample_grouped = df_simulation.groupby(["dataset", "sample"])
    rows = []
    for name, df_g in by_sample_grouped:
        group_rows = [
            name + (cat, df_g.loc[df_g.ncbi.isin([member.rsplit(".",1)[0] for member in cat_members]),'expected_average_coverage'].sum())
            for cat, cat_members in catalogue_groupings.items()]
        rows.extend(group_rows)
    df_catalogue_expect = pd.DataFrame(rows, columns = ["dataset", "sample","catalogue_name", "expected_average_coverage"])
    #df_catalogue_expect.query("dataset == '0_005GB'").to_csv(outdir/"sim.tsv", sep="\t")
    # Collecting counts:
    dir_count_matrices = Path(dir_count_matrices)
    catalouge_count_files = [file for file in dir_count_matrices.glob("*.tsv") if not file.stem=="counts_all"]

    df_count = pd.concat(
        pd.read_csv(file, sep="\t", index_col=0)\
            .reset_index()\
            .rename(columns={'index':'kmer'})\
            .melt(id_vars=["kmer"], var_name='dataset_sample', value_name='count_corrected')\
            .assign(catalogue_name = file.stem)
        for file in catalouge_count_files
    )
    df_count[["dataset","sample"]] = df_count['dataset_sample'].str.split(".", expand=True)


    # Merge expected and conut and add error values
    df_count_combined = pd.merge(df_count, df_catalogue_expect, on = ['dataset', 'sample', 'catalogue_name'])
    df_count_combined["RE"] = (df_count_combined['count_corrected'] - df_count_combined['expected_average_coverage']) / df_count_combined['expected_average_coverage']

   
    # collect MAGinator gene sets
    files_mag_gene_sets = Path(dir_mag_flat).glob("*_kmers.csv")
    df_genes_sets = pd.concat(
        pd.read_csv(file)\
            .assign(catalogue_name = file.stem.rsplit("_kmers",1)[0])
        for file in files_mag_gene_sets
    )

    # Generate visualizations
    catalogues = sorted(df_genes_sets.loc[:,'catalogue_name'].drop_duplicates())
    
    total_view = df_count_combined[["catalogue_name", "dataset"]].drop_duplicates().__len__()

    pbar = tqdm(total =total_view, desc="[Generating plots for each dataset and catalogue]", mininterval=30)
    for name_data, df_data in df_count_combined.groupby("dataset"):
        outdir_sub = outdir/name_data
        outdir_sub.mkdir(parents=True, exist_ok=True)
        for catalogue in catalogues:
            df_gene_i = df_genes_sets.loc[df_genes_sets['catalogue_name']==catalogue,:].dropna(axis=1, how="all")
            df_i = df_data.loc[df_data['catalogue_name']==catalogue, ['kmer', 'dataset_sample', 'RE', "count_corrected", "expected_average_coverage"]]
            
            MAG_genesets = list({"init", "best"}.intersection(set(df_gene_i.columns)))
            hist_data = [
                df_i.loc[:,'RE'].values]
            hist_data.extend([
                df_i.loc[df_i['kmer'].isin(df_gene_i[gs]),'RE'].values
                for gs in MAG_genesets
            ])
            # MAG_genesets
            #     df_i.loc[df_i['kmer'].isin(df_gene_i['init']),'RE'].values,
            #     #df_i.loc[df_i['kmer'].isin(df_gene_i['best']),'RE'].values
            # ]

            MRAE = [np.mean(np.abs(x)) for x in hist_data]

            group_labels = [f'all_{len(df_i.kmer.drop_duplicates())}']
            group_labels.extend([f'MAG_{gs}' for gs in MAG_genesets])

            mse_string = "MRAE: "+ " ".join([f"{label}: {value:0.1f}" for label, value in zip(group_labels, MRAE) ])

            # Create distplot with custom bin_size
            try:
                fig = ff.create_distplot(hist_data, group_labels, bin_size=.1,)
            except Exception as err:
                print(f"[{name_data}:{catalogue}] ERROR "+str(err), file=sys.stderr)
                print(f"[{name_data}:{catalogue}] ERROR "+str(err))
                continue
            fig.update_layout(
                title=f"Histogram with KDE curve showing distribution of relative errors for catelogue groupings.<br><sup>Data: {name_data}:{catalogue} - {mse_string}",
                yaxis_title="Probability Density",
                xaxis_title="Relative Error",
                yaxis1=dict(range=[0,5])
                )
            fig.update_traces(opacity=0.75)

            if catalogue == "NZ_CP053893.1.region002":
                fig.update_xaxes(range=[-1,25])
            else:
                fig.update_xaxes(range=[-1,5])

            file_image = outdir_sub / (catalogue+".png")
            fig.write_image(file_image, width=1000, height=800)
            pbar.update(1)
    pbar.close()

    # #do one full:
    # df_i = df_count_combined
    # hist_data = [
    #             df_i.loc[:,'RE'].values,
    #             df_i.loc[df_i['kmer'].isin(df_gene_i['init']),'RE'].values,
    #             df_i.loc[df_i['kmer'].isin(df_gene_i['best']),'RE'].values
    #         ]

    # MRAE = [np.mean(np.abs(x)) for x in hist_data]

    # group_labels = [
    #     f'all_{len(df_i.kmer.drop_duplicates())}',
    #     'MAG_init',
    #     'MAG_best'
    # ]

    # mse_string = "MRAE: "+ " ".join([f"{label}: {value:0.1f}" for label, value in zip(group_labels, MRAE) ])
    # # Create distplot with custom bin_size
    # fig = ff.create_distplot(hist_data, group_labels, bin_size=.1,)
    # fig.update_layout(
    #     title=f"Histogram with KDE curve showing distribution of relative errors across all datapoints<br><sup>{mse_string}",
    #     yaxis_title="Probability Density",
    #     xaxis_title="Relative Error",
    #     yaxis1=dict(range=[0,5])
    #     )
    # fig.update_traces(opacity=0.75)
    # fig.update_xaxes(range=[-1,25])

    # file_image = outdir / "accross_all.png"
    # fig.write_image(file_image, width=1000, height=800)