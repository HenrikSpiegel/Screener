import sys, os
from typing import List
from Bio import SeqIO
import argparse
#from tqdm import tqdm

def gbk_files_to_fasta(file_paths: List[str], output_file = "data/simulated_data/quantification/input_genomes_gbc.fa"):
    out_fh = open(output_file, "w")
    try:
        for fp in file_paths:
            name = os.path.basename(fp).replace(".gbk","")
            records = SeqIO.parse(fp, "genbank")
            for record in records:
                for feat in record.features:
                    if feat.type == "protocluster":
                        assigned = feat.qualifiers["product"]
                        break
                fasta_entry=f">{name} {assigned} {record.description}\n{record.seq}\n"
                out_fh.write(fasta_entry)
    finally:
        out_fh.close()

def gbk_singlefile_to_fasta(gbk_filename, faa_filename):
    input_handle  = open(gbk_filename, "r")
    output_handle = open(faa_filename, "w")
    try:
        for seq_record in SeqIO.parse(input_handle, "genbank") :

            output_handle.write(f">{seq_record.id} {seq_record.description}\n{seq_record.seq}\n")
    finally: 
        output_handle.close()
        input_handle.close()


if __name__ == "__main__":
    import argparse
    ## Front matter - handle input parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("-i",required=True, help="Antismash output dir")
    parser.add_argument("-o", help="output .fa")
    args = parser.parse_args()

    gbk_files = [os.path.join(args.i, x) for x in os.listdir(args.i) if x.endswith(".gbk") and x != "combined.gbk"]

    if args.o:
        fa_filename = args.o
    else:
        fa_filename = os.path.join(args.i, "combined_bgc.fa")
        print("output: -> ", fa_filename, file=sys.stderr)
    gbk_files_to_fasta(gbk_files, fa_filename)