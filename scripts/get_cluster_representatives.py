from pathlib import Path
import pandas as pd
import json
from Bio import SeqIO

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--family-json", type=Path, required=True)
    parser.add_argument("--similarity-matrix", type=Path, required=True)
    parser.add_argument("--antismash-dir", type=Path, required=True)
    parser.add_argument("--outfile", type=Path, required=True)
    args = parser.parse_args()

    # Get cluster defs:
    cluster_json = Path(args.family_json)
    clusters = json.loads(cluster_json.read_bytes())

    # Get similarity matrix
    similarity_matrix = Path(args.similarity_matrix)
    df_sim = pd.read_csv(similarity_matrix, sep="\t", index_col=0)

    # Determine cluster centers/refs ids
    cluster_reps = {
        clust: df_sim.loc[clusters[clust], clusters[clust]].sum(axis=1).idxmax()
        for clust in clusters.keys()
    }

    # Create a multifasta with the cluster reps.
    antismash_dir = Path(args.antismash_dir)
    outfile = Path(args.outfile)
    outfile.parent.mkdir(parents=True, exist_ok=True)

    fasta_entries = []
    for clust, rep in cluster_reps.items():

        file_gbk = antismash_dir / (rep+".gbk")

        #There should only be a single record
        record = next(SeqIO.parse(file_gbk, "genbank"))

        fasta_entry = f">{clust} {rep}\n{record.seq}"
        fasta_entries.append(fasta_entry)
    outfile.write_text(("\n".join(fasta_entries)))