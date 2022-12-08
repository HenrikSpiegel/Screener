from pathlib import Path
import pandas as pd
import numpy as np
import json
import plotly.express as px
from Bio import SeqIO
from Bio.Seq import Seq
import plotly.graph_objects as go
from plotly.subplots import make_subplots

import configparser

from typing import Literal, Sequence

import sys 

from functions import add_genomic_annotation
from tqdm.notebook import tqdm

def collect_record_dict(fp_antismash_json: Path, fp_id_map: Path, dir_pileups:Path):
    """
    returns dict:
    {
    genome: {
        region: {
            start    : int
            end      : int
            products : [str,]
            seq      : str
            coverage : [int,] #GSA CAMISIM output.
        }
    }
    ...
    }
    """
    combined_region_dict = json.loads(fp_antismash_json.read_bytes())
    
    # Antismash records
    record_dict = dict()
    for record in (pbar := tqdm(combined_region_dict["records"], leave=False)):
        pbar.set_description(f"[Loading Antismash records] : {record['name']}")
        record_dict[record["name"]] = {
            f"region{area_i+1:03d}": {
                "start"   : area['start'],
                "end"     : area["end"],
                "products": area['products'],
                "seq"     : record['seq']['data'][area['start']:area["end"]+1]
            }
            for area_i, area in enumerate(record["areas"])}
    
    

    # Map NCBI id to sam-pileup file.
    dict_pileup_map = dict()
    for line in fp_id_map.read_text().strip().split("\n"):
        ref, id_ = line.split("\t")
        dict_pileup_map[Path(id_).stem.replace("_1", ".1")] = dir_pileups/ (ref+".tsv")
    
    for genome in (pbar := tqdm(record_dict.keys()) ):
        pbar.set_description(f"[Appending pileup coverage] : {genome}")
        df_pileup = pd.read_csv(dict_pileup_map[genome], sep="\t", names=["id","position","coverage"])
        for region in record_dict[genome].keys():
            region_start = record_dict[genome][region]["start"]
            region_end   = record_dict[genome][region]["end"]
            record_dict[genome][region]["coverage"] = df_pileup.loc[region_start:region_end, "coverage"].values.tolist()
              
    return record_dict

def kmerise(seq, k:int=21):
    if len(seq) < k:
        raise ValueError(f"k cannot be larger than lenght of sequence -> k={k} > len(seq)={len(seq)}")
    return [seq[i:i+k] for i in range(len(seq)-k+1)]

def canonicalize(mer:str) -> str:
        mer1 = Seq(mer)
        mer2 = mer1.reverse_complement()
        if mer1 < mer2:
            return str(mer1)
        return str(mer2)


def gen_indexed_geneset(genome, region, record_dict:dict, catalogue_groupings:dict, dir_catalogues:Path, dir_mag_flat:Path):
    record = record_dict[genome][region]
    kmer_index = {
        canonicalize(kmer):i
        for i, kmer in enumerate(kmerise(record["seq"]))
    }
    
    # {familyX: [m1, m2],...} -> {m1: familyX, m2: familyX,..}
    catalogue_belongs = {member: group for group, members in catalogue_groupings.items() for member in members}
    catalogue_name = catalogue_belongs[f"{genome}.{region}"]
    
    #Grab catalogue and refined kemrs
    fp_full_catalogue = dir_catalogues / f"{catalogue_name}.catalogue"
    kmers_full_catalogue = [
        record.seq.__str__()
        for record in SeqIO.parse(fp_full_catalogue, 'fasta')
    ]
    
    fp_mag_gene_sets = dir_mag_flat / (f"{catalogue_name}_kmers.csv")
    
    df_kmers = pd.concat([
        pd.read_csv(fp_mag_gene_sets).melt(value_vars=['init', 'best'], var_name='geneset', value_name='kmer'),
        pd.DataFrame({'geneset':f'all_{len(kmers_full_catalogue)}', 'kmer':kmers_full_catalogue})
    ])
    df_kmers["position"] = pd.Series([kmer_index.get(kmer, None) for kmer in df_kmers["kmer"]],dtype=pd.Int64Dtype())
    
    return df_kmers, kmer_index

def generate_catalogue_coverage_fig(
        genome:str, 
        region:str, 
        record_dict: dict, 
        df_indexed_geneset:pd.DataFrame, 
        antismash_dir:Path = None,
        x_order:Literal['position','coverage'] = 'position'
):

    
    record = record_dict[genome][region]
    
    fig = make_subplots(
        rows=2, cols=1,
        shared_xaxes=True,
        row_width=[0.3, 0.7],
        subplot_titles=["Simulated coverage", "Kmer-set location"],
        vertical_spacing=0.06
    )
    
    if x_order=="position":
        x_plot = range(len(record['coverage']))
        y_plot = record['coverage']
    elif x_order == 'coverage':
        x_plot, y_plot = zip(*[(x,y) for x, y in sorted(
                                        zip(range(len(record['coverage'])), record['coverage']), 
                                        key=lambda pair: pair[1])])
    else:
        raise KeyError("x_order can take ['position', 'coverage'] -> "+ x_order)
        
    fig_cover = px.area(x=x_plot, y=y_plot,
        title=f"Coverage plot for {genome}.{region}",
        labels={
         "x": "Position in Cluster",
         "y": "Coverage (Samtools.Pileup)"
        },
        template="plotly_white"
    )
    
    fig.add_trace(
        fig_cover.data[0]
    )
    
    if antismash_dir and (x_order == 'position'):
        
        file_genebank_record = antismash_dir / (f"{genome}.{region}.gbk")
        gbk_record = next(SeqIO.parse(file_genebank_record, "genbank"))
        add_genomic_annotation(fig_cover, gbk_record)

        fig.update_layout(annotations=fig_cover.layout['annotations'], shapes=fig_cover.layout['shapes'])

    for name, df_geneset in df_indexed_geneset.groupby('geneset'):
        fig.add_trace(
            go.Scatter(
                x = df_geneset['position'],
                y = [name for i in range(len(df_geneset))],
                legendgroup = name,
                name=name,
                mode='markers',
                marker=dict(symbol='line-ns-open')
             ),
            row=2, col=1
        )
    
    if x_order == "coverage":
        fig.update_xaxes(type='category', visible=False)
        fig.update_layout(title=f"{genome}.{region}: Coverage and kmer-set location.<br><sup>Ordered by coverage (BGC: {record['products']})")
    else:
        fig.update_layout(title=f"{genome}.{region}: Coverage and kmer-set location.<br><sup>(BGC: {record['products']})")
    return fig
    
# Lets try plot onto the coverage the kmers with high and low relative error.

def error_correction(values: Sequence[float], k:int=21, error_rate:float=0.03):
    correction_value = (1-error_rate)**k
    corrected_value = [x/correction_value for x in values]
    return corrected_value


def edgeloss_correction(values: Sequence[float], k:int=21, L:int=150):
    correction_value = 1 - ((k-1)/L)
    corrected_value = [x/correction_value for x in values]
    return corrected_value

def combined_correction(values, k, L, error_rate):
    err_corrected = error_correction(values=values, k=k, error_rate=error_rate)

    return edgeloss_correction(
        values=err_corrected,
        k=k, L=L
    )

def get_simulation_overview(file_simulation_overview_fuzzy="data/simulated_data/camisim/*GB/simulation_overview.csv"):
    dfs = []

    fuzz = str(file_simulation_overview_fuzzy)
    if fuzz.startswith("/"):
        summary_glob = Path("/").glob(fuzz[1::])
    else:
        summary_glob = Path().glob(fuzz)
    

    for summary_file in summary_glob:
        df = pd.read_csv(summary_file, sep=",")
        df["dataset"] = summary_file.parent.name
        dfs.append(df)
    df_simulation = pd.concat(dfs)
    
    #Grab config data:
    camisim_config = configparser.ConfigParser()
    #Here we assume only size (ie. readdepth) changes.
    config_file = Path("../data/simulated_data/camisim/configs/").glob("*_config.ini").__next__()
    camisim_config.read(config_file)
    df_simulation["err_type"]=camisim_config.get("ReadSimulator","type")
    try:
        df_simulation["err_profile"]=camisim_config.getfloat("ReadSimulator","profile")
    except:
        df_simulation["err_profile"]=camisim_config.get("ReadSimulator","profile")
    df_simulation["fragment_size"]=camisim_config.getint("ReadSimulator","fragments_size_mean")
    
    return df_simulation

def get_catalogue_errors(fp_simulation_overview, dir_count_matrices, catalogue_groupings):

    df_simulation = get_simulation_overview(fp_simulation_overview)
    df_simulation["readsGB"]  = df_simulation.dataset.str.replace("_",".").str.rstrip("GB").astype(float)

    by_sample_grouped = df_simulation.groupby(["dataset", "sample"])
    rows = []
    for name, df_g in by_sample_grouped:
        group_rows = [
            name + (cat, df_g.loc[df_g.ncbi.isin([member.rsplit(".",1)[0].replace(".","_") for member in cat_members]),'expected_average_coverage'].sum())
            for cat, cat_members in catalogue_groupings.items()]
        rows.extend(group_rows)
    df_catalogue_expect = pd.DataFrame(rows, columns = ["dataset", "sample","catalogue_name", "expected_average_coverage"])

    # Collecting counts:
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
    k = 21 #load from config ...
    L = 148.5
    error_rate = 0.015
    df_count["count_corrected"] = combined_correction(df_count["count"], k, L, error_rate)


    # Merge expected and conut and add error values
    df_count_combined = pd.merge(df_count, df_catalogue_expect, on = ['dataset', 'sample', 'catalogue_name'])
    df_count_combined["RE"] = (df_count_combined['count_corrected'] - df_count_combined['expected_average_coverage']) / df_count_combined['expected_average_coverage']
    df_count_combined["RAE"] = df_count_combined["RE"].abs()
    return df_count_combined

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--catalogues", required=True, type=Path)
    parser.add_argument("--family_dump", required=True, type=Path)
    parser.add_argument("--mag-flat",  required=True, type=Path)
    parser.add_argument("--antismash",  required=True, type=Path)
    parser.add_argument("--simulation-overview", required=True, type=Path)
    parser.add_argument("--count-matrices", required=True, type=Path)
    parser.add_argument("--camisim-id-map", required=True, type=Path)
    parser.add_argument("--pileup-dir", required=True, type=Path)
    parser.add_argument("-o", required=True, type=Path)

    args = parser.parse_args()

    outdir = Path(args.o)
    outdir.mkdir(parents=True, exist_ok=True)

    dir_catalogues  = args.catalogues       #Path("../data/simulated_data/catalogues/catalogues/")
    dir_mag_flat    = args.mag_flat         #Path("../data/simulated_data/MAGinator/screened_flat/")
    dir_antismash   = args.antismash        #Path("../data/simulated_data/antismash/input_genomes/")
    fp_catalogue_groups =  args.family_dump                 #Path("../data/simulated_data/catalogues/family_dump.json")
    catalogue_groupings = json.loads(fp_catalogue_groups.read_bytes())
    catalogue_belongs = {member: group for group, members in catalogue_groupings.items() for member in members}

    # Getting expected values:
    fp_simulation_oveview = args.simulation_overview #"../data/simulated_data/camisim/0_5GB/simulation_overview.csv" #TODO: This could be a place for looping <--
    dir_count_matrices = args.count_matrices ##Path("../data/simulated_data/kmer_quantification/count_matrices/")

    df_catalogue_errors = get_catalogue_errors(fp_simulation_oveview, dir_count_matrices, catalogue_groupings )
    df_average_expected = df_catalogue_errors.groupby("catalogue_name")["expected_average_coverage"].agg("mean").reset_index()
    # NOTE: We are running only on one sample here!
    #df_count_combined
    df_catalogue_errors_agg = df_catalogue_errors.query("dataset=='0_5GB' & sample=='sample_0'").groupby(['catalogue_name','kmer'])["RAE"].agg(["mean","std"]).add_suffix("_RAE").reset_index()

    fp_antismash_json = dir_antismash/ "combined.json"

    fp_id_map = args.camisim_id_map     #Path("../data/simulated_data/camisim/id_map.tsv")

    dir_pileups = args.pileup_dir # Path("../experiments/visualize_catalogue/pileups/") #TODO: This could be a place for looping <--
    record_dict = collect_record_dict(fp_antismash_json, fp_id_map, dir_pileups)

    regions = [(g, r) for g in record_dict.keys() for r in record_dict[g].keys()]
    for genome, region in (pbar := tqdm(regions,mininterval=15)):

        pbar.set_description(f"[{genome}.{region}] Generating data")
        catalogue_name = catalogue_belongs[f"{genome}.{region}"]

        # Generate geneset for kmers belonging to initial 5000 set and the MAG_init and MAG_best
        df_geneset_indexed, kmer_index = gen_indexed_geneset(genome, region, record_dict, catalogue_groupings, dir_catalogues, dir_mag_flat)

        # Extend the geneset with kmers from the initial 5000 annotated to how well they fit the expected value (over all samples)
        df_catalogue_errors_i = df_catalogue_errors_agg.query(f"catalogue_name == '{catalogue_name}'").copy()
        df_catalogue_errors_i['position'] = [kmer_index.get(kmer, None) for kmer in df_catalogue_errors_i["kmer"]]

        geneset_prefix = "Error group: "
        df_catalogue_errors_i["geneset"] = pd.qcut(
            df_catalogue_errors_i.mean_RAE, 
            (0, .25, .5, .75, 1), 
            labels=[geneset_prefix+"<q1", geneset_prefix+"q1-q2", geneset_prefix+">q2", geneset_prefix+">q3"])

        df_geneset_indexed_expanded = pd.concat([
            df_geneset_indexed,
            df_catalogue_errors_i[['kmer','position', 'geneset']]
        ])
        pbar.set_description(f"[{genome}.{region}] Generating plots")
        for fig_type in ["position", "coverage"]:
            fig = generate_catalogue_coverage_fig(
                    genome = genome, 
                    region = region, 
                    record_dict = record_dict, 
                    df_indexed_geneset = df_geneset_indexed_expanded, 
                    antismash_dir=dir_antismash,
                    x_order = fig_type
            )
            fig.update_layout(height=600)
            fp_fig = fp_fig_position = outdir / f"{genome}.{region}_{fig_type}.png"
            fig.write_image(fp_fig, height=600, width=1000)