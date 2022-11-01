from pathlib import Path
import pandas as pd
from Bio import SeqIO
from typing import List

def summarise_blast_output(df_raw):
    df = df_raw.copy()
    df["weighted_pident"] = df.pident * df.length

    summarised = []
    for name, group in df.groupby(["qaccver", "saccver"]):
        qcov = group.qcovs.iloc[0]
        summarised_pident = group.weighted_pident.sum() / group.length.sum()
        summarised.append(name + (qcov,  summarised_pident))
    df_summarised = pd.DataFrame(
        data=summarised,
        columns = ("query", "subject", "query_coverage", "summarised_pident")
    )
    return df_summarised

def similarity_metric(row):
    aligned_bases = row['query_coverage']*row['qlen']*(row['summarised_pident']/100)
    maximum_aligned_bases =   row[['qlen', 'slen']].max()
    return (aligned_bases, aligned_bases / maximum_aligned_bases)

def get_sequence_lengths(mfasta):
    seq_len = {}
    for record in SeqIO.parse(mfasta, "fasta"):
        seq_len[record.name] = len(record.seq)
    return seq_len


if __name__ == "__main__":
    import argparse
    import argparse
    ## Front matter - handle input parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("--blast", required=True, type=Path,  help="blast results outfile")
    parser.add_argument("--fasta", required=True, type=Path, help="query fasta file for blast - required for sequence lengths")
    parser.add_argument("-o", help="output dir - default to directory of blast file.")

    args = parser.parse_args()


    outdir = Path(args.o or args.blast.parent)
 
    out_tsv_summarised          = outdir / "summarised_blast.tsv"
    out_tsv_table_non_symmetric = outdir / "pairwise_table_non_symmetric.tsv"
    out_tsv_table_symmetric     = outdir / "pairwise_table_symmetric.tsv"

    df = pd.read_csv(args.blast, sep="\t")
    df_summarised = summarise_blast_output(df)

    seq_len = get_sequence_lengths(args.fasta)

    df_summarised["qlen"] = [seq_len[x] for x in df_summarised["query"]]
    df_summarised["slen"] = [seq_len[x] for x in df_summarised["subject"]]
    
    df_summarised[["aligned_bases", "similarity"]] = df_summarised.apply(similarity_metric, axis=1, result_type="expand")
    df_summarised.to_csv(out_tsv_summarised, sep="\t", index=False)

    df_wide = df_summarised.pivot_table(index="query", columns="subject", values="similarity")
    df_wide.to_csv(out_tsv_table_non_symmetric, sep="\t", index=True)

    values = df_wide.values
    sym_values = values = (values + values.T) / 2
    df_wide_sym = pd.DataFrame(sym_values, columns=df_wide.columns, index=df_wide.index)
    df_wide_sym.to_csv(out_tsv_table_symmetric, sep="\t", index=True)