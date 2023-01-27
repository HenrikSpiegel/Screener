"""
Extract and combined the .gbk files from multiple antismash runs.

Generates multiple temporary combine files (for each run) and concats them).

"""

import argparse
from pathlib import Path
from multiprocessing import Pool
import subprocess
from typing import List
import shutil
from Bio import SeqIO


def parse_dir(directory: Path):
    call = [
        "python", "scripts/antismash_as_fasta.py",
        "-i", directory
    ]
    subprocess.run(call, check=True)
    return 0

# def cleanup_names(fasta_file:Path):
#     outfile = fasta_file.with_name(fasta_file.name+".corrected")
#     sample_name = fasta_file.name.split(".")[0].split("_")[1]
#     with open(outfile, 'w') as corrected:
#         for record in SeqIO.parse(fasta_file, "fasta"): 
#             contig = record.id.split("_")[0]
#             region = record.id.split(".")[-1]
#             new_id = sample_name+"_"+contig+"."+region
#             new_desc = new_id + " " +record.description.split(" ", 1)[-1]
#             record.id = new_id
#             record.description = new_desc
#             SeqIO.write(record, corrected, 'fasta')
#     return outfile

def copy_and_rename(antismash_dir: Path, outdir = Path):
    sample_name = antismash_dir.name.split(".")[0].split("_")[1]
    for fp_gbk in antismash_dir.glob("*region*.gbk"):
        for record in SeqIO.parse(fp_gbk, "genbank"):
            contig = fp_gbk.stem.split("_")[0]
            region = fp_gbk.stem.split(".")[-1]
            new_id = sample_name+"_"+contig+"."+region
            new_desc = new_id + " " +record.description.split(" ", 1)[-1]
            record.id = new_id
            #record.name = new_id
            record.description = new_desc

            outfile = outdir / (new_id+".gbk")
            SeqIO.write(record, outfile, 'genbank')

def collect(temp_files: List[Path], outfile: Path):
    with open(outfile,'wb') as wfd:
        for f in temp_files:
            with open(f,'rb') as fd:
                shutil.copyfileobj(fd, wfd)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--antismash-dirs", required=True, type=Path, nargs="+")
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    antismash_dirs = args.antismash_dirs

    parsing_args = [(as_dir, outdir) for as_dir in antismash_dirs]    

    print(f"Running parsers for ({len(antismash_dirs)}) dirs")
    with Pool(args.threads) as p:
        p.starmap(copy_and_rename, parsing_args)

    print("Collecting files")
    parse_dir(outdir)

    print(f"finished -> {outdir}")






