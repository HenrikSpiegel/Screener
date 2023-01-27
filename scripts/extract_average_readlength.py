from pathlib import Path
import gzip
import sys
from Bio import SeqIO
from typing import List
import json

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
    try:
        avg_readlength = combined_length / n_reads
    except ZeroDivisionError:
        print(files, file=sys.stderr)
        raise
    return avg_readlength

if __name__ == "__main__":
    import argparse
    from multiprocessing import Pool
    

    parser = argparse.ArgumentParser()
    parser.add_argument("--top-dir", required=True, type=Path, help="Top dir containing the fuzzy dataset paths")
    parser.add_argument("--fuzzy-dataset-names", default="*GB/sample*", type=str, help="Fuzzy path for dataset subdirs ['*GB/sample*']")
    parser.add_argument("-o", type=Path, help="redirect out from --top-dir/average_lenghts.json")
    parser.add_argument("--threads", type=int, default=5, help="Maximum number of workers started to count individual datasets.[5]")
    args = parser.parse_args()


    dir_preprocessed   = args.top_dir
    fuzzy_dataset_name = args.fuzzy_dataset_names

    out_json = args.o or dir_preprocessed/"average_lengths.json"

    #fuzzy_readname = "*reads.fq*"

    dataset_paths = list(dir_preprocessed.glob(fuzzy_dataset_name))
    if not dataset_paths:
        raise FileExistsError(f"{dir_preprocessed}/{fuzzy_dataset_name}")

    dataset_read_files = [list(ds.iterdir()) for ds in dataset_paths]

    print(f"Checking read files using: ({args.threads}) readers", file=sys.stderr)
    with Pool(args.threads) as p:
        avg_lengths = p.map(get_average_length, dataset_read_files)

    dict_out = {
        ds.as_posix():avg_len
        for ds, avg_len in zip(dataset_paths, avg_lengths)
    }

    print(f"output -> {out_json}", file=sys.stderr)
    out_json.write_text(json.dumps(dict_out))
    