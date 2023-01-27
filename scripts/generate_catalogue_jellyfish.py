import json
from Bio import SeqIO
from pathlib import Path
import os
import subprocess
from multiprocessing import Pool
import pandas as pd
from tqdm import tqdm


def jelly_count(fasta: Path, out_name: Path) -> None:
    cmd = f"jellyfish count -m 21 -s 5G -t 4 -C -o {out_name} {fasta}"
    subprocess.run(cmd.split(" "), check=True)


def jelly_dump(jelly_db: Path, outfile: Path):
    cmd = f"jellyfish dump -o {outfile} -L 1 {jelly_db}"
    subprocess.run(cmd.split(" "), check=True)

def jelly_query(jelly_db: Path, fasta_file: Path, outfile: Path):
    cmd = f"jellyfish query -s {fasta_file} -o {outfile} {jelly_db}"
    subprocess.run(cmd.split(" "), check=True)


def gbk2fasta_combine(gbk_files, outfasta):
    with open(outfasta, "w") as output_handle:
        for gbk_file in gbk_files:
            with open(gbk_file, "r") as input_handle:
                gbk_records = SeqIO.parse(input_handle, "genbank")
                SeqIO.write(gbk_records, output_handle, "fasta")
                
def find_distinct_kmers(dump_file: Path, query_file: Path):
    keep_records = []
    dump_records = SeqIO.parse(dump_file, "fasta")

    query_fh = open(query_file, "r")
    try:
        for dump_record, query_line in zip(dump_records, query_fh):
            q_seq, q_count = query_line.strip().split(" ") 
            if not str(dump_record.seq) == q_seq:
                raise RuntimeError(f"Missmatch: {q_seq} != {dump_record.seq}")
            if q_count == dump_record.id:
                keep_records.append(dump_record)
    except Exception:
        raise
    finally:
        query_fh.close()
    
    keep_records_sorted = sorted(keep_records, key= lambda seq: int(seq.id), reverse=True)
    return keep_records_sorted

def generate_meta_file(sorted_distinct, keep_kmers:int, cluster_size:int, outfile: Path) -> None:
    df_out = pd.DataFrame(
        {
            'kmer': str(r.seq),
            'freq': int(r.id)/cluster_size
        }
        for r in sorted_distinct[0:keep_kmers]  
    ).assign(size=21)
    
    df_out.to_csv(outfile, sep="\t", index=False, header=False)
    
def generate_catalogue_fasta(sorted_distinct, keep_kmers:int, cluster_size:int, outfile: Path):
    
    with open(outfile, "w") as output_handle:
        for i, record in enumerate(sorted_distinct):
            new_id = f"kmer.{i} prevalence {int(record.id)/cluster_size:.2f}"
            record.id = new_id
            SeqIO.write(record, output_handle, "fasta")
            if i == keep_kmers-1:
                break

def run_cluster(wd, cluster_name, cluster_members, keep_lenght) -> None:
    
    dir_jelly = wd/"jelly"
    dir_jelly.mkdir(parents=True, exist_ok=True)

    dir_fasta = wd/"fasta"
    dir_fasta.mkdir(parents=True, exist_ok=True)

    dir_catalogue = wd / "catalogues"
    dir_meta = wd / "metafiles"

    cluster_fa      = dir_fasta / (cluster_name+".fa")
    cluster_jelly_db      = dir_jelly / (cluster_name+".jf")
    cluster_jelly_dump    = dir_jelly / (cluster_name+".dump.fa")
    cluster_jelly_qcount  = dir_jelly / (cluster_name+".query_all.csv")

    cluster_meta = dir_meta / (cluster_name+".meta")
    cluster_catalogue = dir_catalogue / (cluster_name+".catalogue")

    # Generate multifasta for cluster
    gbk_files = [dir_antismash / (member+".gbk") for member in cluster_members]
    gbk2fasta_combine(gbk_files, cluster_fa)

    jelly_count(cluster_fa, cluster_jelly_db)
    jelly_dump(cluster_jelly_db, cluster_jelly_dump)

    # Compare cluster dump with counts of all
    jelly_query(fp_count_all, cluster_jelly_dump, cluster_jelly_qcount)
    sorted_distinct_kmers = find_distinct_kmers(cluster_jelly_dump, cluster_jelly_qcount)

    # Generate catalogue fasta and metafiles.
    generate_meta_file(
        sorted_distinct = sorted_distinct_kmers,
        keep_kmers = keep_lenght,
        cluster_size = len(cluster_members), 
        outfile = cluster_meta)

    generate_catalogue_fasta(
        sorted_distinct_kmers, 
        keep_kmers=keep_lenght, 
        cluster_size=len(cluster_members), 
        outfile=cluster_catalogue)

if __name__ == "__main__":
    import argparse

    desc = """\
Generate distinct sets of kmers based on a clusters as defined in a json'ed dictionary:
{cluster_name: [mem1, mem2, ....], ...}

Assumes jellyfish is availble from path.
"""

    parser= argparse.ArgumentParser(description=desc)
    parser.add_argument("--cluster-def-json", type=Path, required=True, help="Path to .json describing cluster membership")
    parser.add_argument("--dir-antismash", type=Path, required=True, help="path to dir containing the antismash output - .gbk files must match member names in cluster-def-json ")
    parser.add_argument("--wd", type=Path, required=True, help="Working dir for script / jellyfish dumps etc.")
    parser.add_argument("--n-kmers", type=int, default=10000, help="Maximum number of distinct kmers per cluster [10000]")
    args = parser.parse_args()

    fp_clusters = args.cluster_def_json     # Path("data/ibdmdb/mcl_clustering/out.blast_result.mci.I40.json")
    clusters = json.loads(fp_clusters.read_bytes())

    dir_antismash = args.dir_antismash                #  Path("data/ibdmdb/antismash")
    fasta_all = dir_antismash / "combined_bgc.fa"

    wd = args.wd        #Path("data/ibdmdb/catalogues/")
    dir_jelly = wd /"jelly"

    keep_lenght = args.n_kmers

    # Generate count of global kmers:
    fp_count_all = dir_jelly/"all.jf"
    if not fp_count_all.is_file():
        jelly_count(fasta=fasta_all, out_name=fp_count_all)

    # Note this could easily be a pool of subjobs instead of a loop to make it parallele

    args = [
        (wd, cluster_name, cluster_members, keep_lenght)
        for cluster_name, cluster_members in clusters.items()
        ]

    with Pool(5) as p:
        p.starmap(run_cluster, tqdm(args, total=len(args)), chunksize=2)
    #for arg in args:
    #    run_cluster(*arg)
    #    break


