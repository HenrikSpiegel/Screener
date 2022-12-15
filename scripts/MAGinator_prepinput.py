from pathlib import Path
import numpy as np
import pandas as pd
import sys
import shutil


# geneID_collectionID.tsv
def generate_MAGinator_files(dir_catalogue, dir_MAGinator):
    dir_MAGinator = Path(dir_MAGinator)
 
    file_geneID_collectionID = dir_MAGinator / 'clusters/gene_lists/geneID_collectionID.tsv'
    file_fasta_fai = dir_MAGinator / 'genes/all_genes_nonredundant.fasta.fai'
    file_clustername_order = dir_MAGinator /"collectionID_order.txt"

    file_geneID_collectionID.parent.mkdir(parents=True, exist_ok=True)
    file_fasta_fai.parent.mkdir(parents=True, exist_ok=True)

    file_geneID_collectionID.unlink(missing_ok=True)
    file_fasta_fai.unlink(missing_ok=True)

    clustername_order = []
    fh_geneID_collectionID = file_geneID_collectionID.open(mode="a")
    fh_fasta_fai = file_fasta_fai.open(mode="a")
    for i, catalogue in enumerate(dir_catalogue.iterdir()):
        clustername_order.append(catalogue.stem)
        kmers = [x.strip() for x in catalogue.read_text().split("\n") if not x.startswith(">")]

        fh_geneID_collectionID.write(
            "\n".join(f"{x}\t{i}" for x in kmers) + '\n' 
            
        )
        fh_fasta_fai.write(
            "\n".join(f"{x}\t1" for x in kmers) + '\n'
            
        )
    fh_geneID_collectionID.close()
    fh_fasta_fai.close()
                                     
    file_clustername_order.write_text("\n".join(clustername_order))
    print(f"Created geneID_collectionID with ({i}) collections -> {file_geneID_collectionID}", file=sys.stderr)
    print(f"Created all_genes_nonredundant.fasta.fai -> {file_fasta_fai}", file=sys.stderr)

    #Add required but unused files.
    unused_dir = dir_MAGinator / 'unused'
    unused_dir.mkdir(parents=True, exist_ok=True)
    (unused_dir/"contigs.fasta").touch()
    (unused_dir/"clusters.tsv").touch()
    read_file = unused_dir/"fake_reads.fa"
    read_file.write_text("fakesample1, fakesample2")

    tab_file = (dir_MAGinator/'tabs/metagenomicspecies.tab')
    tab_file.parent.mkdir(exist_ok=True, parents=True)
    tab_file.touch()

    vamb_file = unused_dir/"fake_vamb.file"
    vamb_file.touch()
    par_file = unused_dir/'parameters.tab'
    par_file.write_text(f"reads {read_file.resolve()}\nvamb_clusters {vamb_file.resolve()}")



    
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--catalogues", required=True, type=Path, help="directory containing .catalogue files")
    parser.add_argument("--count-matrix", required=True, type=Path, help="filepath for countmatrix for all catalogues and all samples.")
    parser.add_argument("--min-dataset", type=float, help="If set will only use datasets >= min-dataset. (0_05GB.samplex -> 0.05gb)")
    parser.add_argument("-o", required=True, type=Path, help="top directory for MAGinator data workingdir")
    args = parser.parse_args()

    dir_catalogue = args.catalogues # Path("../data/simulated_data/catalogues/catalogues/")
    file_counts_all = args.count_matrix #Path("../data/simulated_data/kmer_quantification/count_matrices/counts_all.tsv")

    dir_MAGinator = args.o # Path("../data/simulated_data/MAGinator")
    file_MAGinator_count = dir_MAGinator/"genes/small_gene_count_matrix.tsv"
    file_MAGinator_count.parent.mkdir(parents=True, exist_ok=True)

    generate_MAGinator_files(dir_catalogue, dir_MAGinator)
    if not args.min_dataset:
        shutil.copyfile(file_counts_all, file_MAGinator_count)
    else:
        df_counts_all = pd.read_csv(file_counts_all, index_col=0)
        col_ds_sizes = np.array([float(x.split(".")[0].replace("GB","").replace("_",".")) for x in df_counts_all.columns])
        cols_include = df_counts_all.columns[col_ds_sizes>args.min_dataset]
        df_counts_include = df_counts_all.loc[:,cols_include]
        df_counts_include.to_csv(file_MAGinator_count, sep="\t")
        
    print(f"Created countmatrix at -> {file_MAGinator_count}", file=sys.stderr)

