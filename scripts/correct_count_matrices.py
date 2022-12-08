from pathlib import Path
import pandas as pd

from Bio import SeqIO

from typing import List, Sequence
import gzip
from tqdm import tqdm

def get_average_length(files: List[Path]):
    combined_length = 0
    n_reads = 0
   
    for file in files:
        if not file.exists():
            raise IOError(f"{file} not found")
        with gzip.open(file, "rt") as handle:
            for record in SeqIO.parse(handle, "fastq"):
                combined_length += record.seq.__len__()
                n_reads += 1
    avg_readlength = combined_length / n_reads
    return avg_readlength

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

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--counts", required=True, type=Path, help="dir for raw count matrices")
    parser.add_argument("--reads-dir", required=True, type=Path, help="Dir holding processed reads (or raw if not processing)")
    parser.add_argument("--fuzzy-dataset-names", type=str, default='*GB/sample_*', help="fuzzy names for datasets in reads_dir, assumes parent*/sample* style naming.['*GB/sample_*']")
    parser.add_argument("-k", default=21, type=int)
    parser.add_argument("--est-read-err", default=0.015, type=float, help="Estimated average per base read error.")
    
    parser.add_argument("-o", required=True, type=Path, help="outdir for corrected matrices.")
    args=parser.parse_args()

    ##

    dir_counts_raw = args.counts        # WD/"kmer_quantification/count_matrices"
    dir_counts_corrected = args.o           #WD/"kmer_quantification/count_matrices_corrected"

    dir_preprocessed = args.reads_dir #WD/"preprocessed"
    fuzzy_dataset_name = args.fuzzy_dataset_names ##"*GB/sample*"

    kmer=args.k
    est_error=args.est_read_err

    dir_counts_corrected.mkdir(parents=True, exist_ok=True)

    ##
    #print(dir_preprocessed)
    #print(list(dir_preprocessed.glob(fuzzy_dataset_name)))

    # Get average readlength: Quite expensive - could perhaps be extracted from the preprocess stderr/stdout???
    sample_avg_readlength = {}
    datasets = dir_preprocessed.glob(fuzzy_dataset_name)
    if not datasets:
        raise FileNotFoundError(dir_preprocessed/fuzzy_dataset_name)
    for dataset in datasets:      #(pbar := tqdm(list((dir_preprocessed).glob(fuzzy_dataset_name)))):
        sample_short = dataset.parent.name+"."+dataset.name
        #pbar.set_description(f"[{sample_short}] Getting average processed readlength")
        
        files = dataset.glob("*reads.fq.gz")
        avg_read = get_average_length(files)
        sample_avg_readlength[sample_short] = avg_read

    for fp_raw in dir_counts_raw.glob("*.tsv"):
        filename = fp_raw.name
        df_count = pd.read_csv(fp_raw, sep="\t", index_col=0)
        
        df_count_corrected = pd.DataFrame(
            index = df_count.index,
            data = {
                column: combined_correction(df_count[column], k=kmer, error_rate=est_error, L=sample_avg_readlength[column])
                for column in df_count.columns
            }
        )
        outfile = dir_counts_corrected/filename
        df_count_corrected.to_csv(outfile, sep="\t", index=True)