#!/usr/bin/env python3

## Python script for pulling files from ncbi.

try:
    import requests
    import os, sys
    from typing import List
    import shutil
    import pathlib
except:
    print("Ensure modules are loaded: suggested (tools, module load anaconda3/2020.07", file=sys.stderr)


def fetch_fasta_from_ncbi(ncbi_id, outdir = None, overwrite: bool=False, verbose: bool=False) -> str:
    
    
    entrez_api_url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={ncbi_id}&rettype=fasta&retmode=text'
    local_filename = ncbi_id.replace(".","_")+".fa"
    
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

def fetch_fastas_async(ncbi_ids: List[str], outdir: str=None, overwrite: bool=False, verbose: bool=False) -> List[str]:
    loaded_files = []
    for ncbi_id in ncbi_ids: #todo: run as async group
        filepath = fetch_fasta_from_ncbi(ncbi_id = ncbi_id, outdir = outdir, overwrite=overwrite, verbose=verbose)
        loaded_files.append(filepath)
    return loaded_files

if __name__ == "__main__":
    # Is run as script:
    
    import argparse, sys
    ## Front matter - handle input parsing
    parser = argparse.ArgumentParser()
    parser.add_argument('--ncbi_ids', default=sys.stdin, nargs='+', help="list of ids to be fetched, can be sequence stdin stream")
    parser.add_argument('--outdir', default=None, help="target for outdir")
    parser.add_argument("--overwrite", action="store_true", help="overwrite existing .fasta files")
    parser.add_argument("--verbose", action="store_true", help="Increased logging to stderr")

    args = parser.parse_args()
    #print(args.ncbi_ids.name=="<stdin>")
    #print([entry.strip() for line in args.ncbi_ids.readlines() for entry in line.split() ])

    # get ids:
    if type(args.ncbi_ids) == list:
        ncbi_ids = args.ncbi_ids
        # TODO: Add some sort of check?
    elif args.ncbi_ids.name=="<stdin>":
        try:
            ncbi_ids = [entry.strip() for line in args.ncbi_ids.readlines() for entry in line.split()]
        except:
            print("Failed to load entries from file. Please arrange as 1 entry per line.", file=sys.stderr)
            raise
    else:
        raise KeyError("Failed reading __ncbi_ids")
        
    ## Main matter - run the loading.
    loaded_files = fetch_fastas_async(ncbi_ids = ncbi_ids,
                                      outdir = args.outdir,
                                      overwrite = args.overwrite,
                                      verbose=args.overwrite)
    ## combine all the files.
    combined_fp = os.path.join(os.path.dirname(loaded_files[0]), "combined.fa")
    with open(combined_fp,'wb') as wfd:
        for f in loaded_files:
            with open(f,'rb') as fd:
                shutil.copyfileobj(fd, wfd)

    print(*loaded_files, sep="\n")
