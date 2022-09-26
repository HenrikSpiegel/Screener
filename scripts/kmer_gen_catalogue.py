#!/usr/bin/env python3
import pathlib
import sys
import os
from Bio.Seq import Seq

def kmerise(seq, k:int=21):
    if len(seq) < k:
        raise ValueError(f"k cannot be larger than lenght of sequence -> k={k} > len(seq)={len(seq)}")
    return [seq[i:i+k] for i in range(len(seq)-k+1)]

def canonicalize(mer:str) -> str:
    mer1 = Seq(mer)
    mer2 = mer1.reverse_complement()
    if mer1 < mer2:
        return str(mer1)
    return str(mer2)

def generate_raw_kmer_catalogue(fasta: str, outdir: str="", k:int=21) -> int:
    catalogues = 0
    header = None
    seq    = ""
    for line in open(fasta, "r"):
        if line.startswith(">"):
            #Write buffer
            if header:
                filename = header[1::].split(" ",1)[0] +".catalogue_raw.fa" #Get first word in 
                fp_out = os.path.join(outdir, filename)
                print("writting -> "+fp_out, file=sys.stderr)
                with open(fp_out, "w") as fh:
                    kmers_cannonical = {canonicalize(x) for x in kmerise(seq, k)}
                    entries = [f">kmer{i}\n{kmer}" for i, kmer in enumerate(kmers_cannonical)]
                    fh.write("\n".join(entries))
                    catalogues += 1
                
            #reset buffer
            header = line.strip()
            seq = ""
        else:
            seq += line.strip()
    #Write final buffer
    if header:
        filename = header[1::].split(" ",1)[0] +".catalogue_raw.fa" #Get first word in header
        fp_out = os.path.join(outdir, filename)
        print("writting -> "+fp_out, file=sys.stderr)
        with open(fp_out, "w") as fh:
            kmers_cannonical = {canonicalize(x) for x in kmerise(seq, k)}
            entries = [f">kmer{i}\n{kmer}" for i, kmer in enumerate(kmers_cannonical)]
            fh.write("\n".join(entries))
            catalogues += 1
    return catalogues

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('-o', "--outdir", help='output directory', required=True)
    parser.add_argument('-f', "--fastas", nargs='+', help='path(s) to fastafiles to catalogue (multi fasta is split)', required=True)
    parser.add_argument('-k', default = 15, type=int, help="lenght of mers.")
    args = parser.parse_args()

    pathlib.Path(args.outdir).mkdir(parents=True, exist_ok=True)
    sys.stderr.write(f"Generating catalogue from ({len(args.fastas)}) fastafiles. \n")
    i = 0
    for fasta in args.fastas:
        i += generate_raw_kmer_catalogue(fasta=fasta, outdir=args.outdir, k=args.k)
        
    sys.stderr.write(f"Finished generating ({i}) catalogues\n")
