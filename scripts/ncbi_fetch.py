#!/usr/bin/env python3

## Python script for pulling files from ncbi.

try:
    import requests
    import os, sys
    from typing import List
    import shutil
    import pathlib
    from pathlib import Path
    import pandas as pd
except:
    print("Ensure modules are loaded: suggested (tools, module load anaconda3/2020.07", file=sys.stderr)


def fetch_fasta_from_ncbi(ncbi_id, outdir = None, overwrite: bool=False, verbose: bool=False) -> str:
    
    
    entrez_api_url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={ncbi_id}&rettype=fasta&retmode=text'
    local_filename = ncbi_id+".fa"#.replace(".","_")+".fa"
    
    if outdir:
        pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)
        local_filehandle = os.path.join(outdir, local_filename)
    else:
        local_filehandle = local_filename
   
    #Check if file already exists:
    if not overwrite and os.path.isfile(local_filehandle):
        if verbose:
            print(f"{ncbi_id} already in outdir -> Skipping", file=sys.stderr)
        return os.path.abspath(local_filehandle)
        
    
    # NOTE the stream=True parameter below
    print("Pulling file with id: " + ncbi_id, file=sys.stderr)
    with requests.get(entrez_api_url, stream=True) as r:
        r.raise_for_status()
        with open(local_filehandle, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192): 
                # If you have chunk encoded response uncomment if
                # and set chunk_size parameter to None.
                #if chunk: 
                f.write(chunk)
    return os.path.abspath(local_filehandle)

def fetch_fastas(ncbi_ids: List[str], outdir: str=None, overwrite: bool=False, verbose: bool=False) -> List[str]:
    loaded_files = []
    for ncbi_id in ncbi_ids: 
        filepath = fetch_fasta_from_ncbi(ncbi_id = ncbi_id, outdir = outdir, overwrite=overwrite, verbose=verbose)
        loaded_files.append(filepath)
    return loaded_files

if __name__ == "__main__":
    # Is run as script:
    
    import argparse, sys
    ## Front matter - handle input parsing
    parser = argparse.ArgumentParser()
    parser.add_argument('--acc_id', nargs='+', type=str, help="list of ids to be fetched")
    parser.add_argument('--tax_id', nargs='+', type=int, help="list of tax-ids matching acc_id order")
    parser.add_argument('--outdir', required=True, type=Path, help="target for outdir")
    parser.add_argument("--overwrite", action="store_true", help="overwrite existing .fasta files")
    parser.add_argument("--verbose", action="store_true", help="Increased logging to stderr")

    args = parser.parse_args()


    # get ids:
    if type(args.acc_id) == list:
        acc_id = args.acc_id
    else:
        acc_id = [args.acc_id]

    
        
    ## Main matter - run the loading.
    loaded_files = fetch_fastas(ncbi_ids = acc_id,
                                      outdir = args.outdir,
                                      overwrite = args.overwrite,
                                      verbose=args.overwrite)
    ## combine all the files.
    combined_fp = os.path.join(os.path.dirname(loaded_files[0]), "combined.fa")
    with open(combined_fp,'wb') as wfd:
        for f in loaded_files:
            with open(f,'rb') as fd:
                shutil.copyfileobj(fd, wfd)

    acc_to_tax = {a:t for a,t in zip(args.acc_id, args.tax_id)}

    ## write _selected.tsv overview.
    df_selected = pd.DataFrame([
        {
            #"genus": " ".join(Path(file).read_text().split("\n")[0][1::].split(" ")[1:3]),
            "assembly_accession": Path(file).stem,
            "taxid": acc_to_tax[Path(file).stem]
        }
        for file in loaded_files
    ])
    df_selected.to_csv(args.outdir/"ncbi_selected.tsv",sep="\t",index=False)

    print(*loaded_files, sep="\n")
