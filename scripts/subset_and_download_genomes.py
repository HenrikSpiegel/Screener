import gzip
from io import BytesIO
from pathlib import Path
import pandas as pd
import numpy as np
from tqdm import tqdm
import requests
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--genera-tables", type=Path, nargs="+", help="Space seperated filepaths")
parser.add_argument("-o", type=Path, required=True)
args = parser.parse_args()

outdir      = Path(args.o)
genomes_dir = outdir / "genomes"
genomes_dir.mkdir(parents=True, exist_ok=True)

if args.genera_tables:
    genera_files = [Path(f) for f in args.genera_tables]
else:
   genus_dir = Path("data/simulated_data_300genomes/input_genomes/genera_tables")
   genera_files = [x for x in genus_dir.iterdir() if x.name != "assembly_summary.txt"]

assert(all(file.is_file() for file in genera_files))

for genus_table_file in (pbar := tqdm( genera_files ) ) :
    pbar.set_description(f"[{genus_table_file.stem}] Selecting subset")

    genus_assemblies_all = pd.read_csv(genus_table_file, sep="\t", low_memory=False
    )
    genus_assemblies_complete = genus_assemblies_all.query("assembly_level == 'Complete Genome'")
    
    #Ensure uniqueness and get the latest version if multiple.
    genus_subset = pd.DataFrame(genus_assemblies_complete.copy()["assembly_accession"])
    genus_subset[["assembly","accession"]] = genus_subset.assembly_accession.str.rsplit(".",1, expand=True)
    genus_subset.sort_values("assembly_accession").drop_duplicates("assembly", keep="last", inplace=True)
    genus_assemblies_unique = genus_subset["assembly_accession"].tolist()
    np.random.seed(1337)

    # random subset if there is more than 100
    n_selected = min(100, len(genus_assemblies_unique))
    selected_assemblies = np.random.choice(genus_assemblies_unique, size=n_selected, replace=False)
    df_selected_assemblies = genus_assemblies_all.loc[genus_assemblies_all.assembly_accession.isin(selected_assemblies),:].reset_index(drop=True)
    
    pbar.set_description(f"[{genus_table_file.stem}] writing subsetted assembly")
    outfile = outdir / (genus_table_file.stem.replace(" ","_") + "_selected.tsv")
    df_selected_assemblies.to_csv(outfile, sep="\t", index=False)
    
    pbar.set_description(f"[{genus_table_file.stem}] Downloading genomes")
    for row in (pbar2 := tqdm(df_selected_assemblies.itertuples(), leave=False) ):
        assembly_acc = row.assembly_accession
        ftp_name =  row.ftp_path.rsplit("/",1)[1]
        pbar2.set_description(f"[{assembly_acc}] downloading....")
        full_ftp_path = row.ftp_path + "/" + ftp_name+"_genomic.fna.gz"
        
        #assembly_short_name = assembly_name.rsplit("_", 1)[0]
        outfile = genomes_dir / (assembly_acc + ".fna.gz")

        bio = BytesIO(requests.get(full_ftp_path).content)
        bio.seek(0)
        #Modify the header so it uses the assembly_ID.version instead of the ncbi_acc
        with gzip.open(bio, "rb") as f:
            data = f.read()
        new_head_begin = f">{assembly_acc} ".encode()
        new_data = new_head_begin + data[1::]
        
        with gzip.open(outfile, "wb") as f:
            f.write(new_data)
