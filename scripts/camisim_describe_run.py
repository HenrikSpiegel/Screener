import glob
import sys
import pandas as pd
from Bio import SeqIO
import os
from pathlib import Path

def generate_simulation_overview(dir_dataset:str):
    dir_dataset = Path(dir_dataset)
    readsGB = float(dir_dataset.name.replace("GB","").replace("_",".")) #very brute - could be nicer
    fp_map = dir_dataset/ "internal/genome_locations.tsv"
        
    cam2genom = {}
    genome_size = {}
    for line in open(fp_map, "r"):
        genomeid, fp = line.strip().split("\t")
        ncbi_id = Path(fp).stem
        record = next(SeqIO.parse(fp, "fasta"))
        cam2genom[genomeid] = ncbi_id
        cam2genom[ncbi_id] = genomeid
        genome_size[ncbi_id] = record.seq.__len__()
        
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
    #parser.add_argument('-f', '--combined_genomes', required=True, help="file containing the combined genomes")
    args = parser.parse_args()

    df = generate_simulation_overview(args.directory)
    out = os.path.join(args.directory, "simulation_overview.csv")
    df.to_csv(out, index=False)
    print("summary_written -> "+out, file=sys.stderr)