import glob
import sys
import pandas as pd
from Bio import SeqIO
import os

def generate_simulation_overview(dir_dataset:str):

    readsGB = float(os.path.basename(dir_dataset.rstrip("/").replace("GB","").replace('_','.'))) #very brute - could be nicer

    #Get genome sizes
    fp_combined_input_genomes = "data/simulated_data/input_genomes/combined.fa" #TODO: fix hardcoded?
    genome_size = {}
    for record in SeqIO.parse(fp_combined_input_genomes, "fasta"):
        genome_size[record.name.replace(".","_")] = len(record.seq)
        
    #Get mapping between camisim id and ncbi ids
    fp_map = os.path.join(dir_dataset, "internal/genome_locations.tsv")
    cam2genom = {}
    for line in open(fp_map, "r"):
        genomeid, fp = line.strip().split("\t")
        ncbi_id = os.path.basename(fp).rsplit(".",1)[0]
        cam2genom[genomeid] = ncbi_id
        cam2genom[ncbi_id] = genomeid
        
    #Run over the samples and calculate the expected coverage from genome size and distribution.
    dfs = []
    for fp_dist in glob.glob(os.path.join(dir_dataset, "distributions/*")):
        df = pd.read_csv(fp_dist, sep="\t", names=["id", "distribution"])
        df["sample"] = "sample_" + fp_dist.rsplit("_",1)[1].replace(".txt","")
        df["ncbi"] = [cam2genom[x] for x in df.id]
        df["genome_size"] = [genome_size[x] for x in df.ncbi]
        df["weighted_dist"] = df.distribution * df.genome_size
        df["readshare"] = (df.weighted_dist/df.weighted_dist.sum())*(readsGB*10**9)
        df["expected_average_coverage"] = df.readshare / df.genome_size

        dfs.append(df)
    df_total = pd.concat(dfs)
    return df_total

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-d','--directory', required=True, help="Directory holding camisim output.")
    args = parser.parse_args()

    df = generate_simulation_overview(args.directory)
    out = os.path.join(args.directory, "simulation_overview.csv")
    df.to_csv(out, index=False)
    print("summary_written -> "+out, file=sys.stderr)