# Get the average readlength of fastqc file(s).
import sys
from typing import List
from Bio import SeqIO
from pathlib import Path
import gzip

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

if __name__ == "__main__":
    desc = "Calculates the average readlength of fastq file(s) and returns as stdout."
    import argparse
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-f', '--files', nargs="+", required=True, help="Readfiles")
    
    args = parser.parse_args()
    file_paths = [Path(f) for f in args.files]
    print(f"Getting average length across ({len(file_paths)} file(s)", file=sys.stderr)
    average_length = get_average_length(file_paths)
    print(average_length)
