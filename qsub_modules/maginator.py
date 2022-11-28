from qsub_modules.base import Base
import os
import pathlib
from pathlib import Path
import glob
import logging
from typing import List, Union

# #modules required: tools anaconda3/2021.05 mamba-org/mamba/0.24.0


# dir_MAGinator="/home/projects/dtu_00009/people/henspi/git/MAGinator/"
# snake_pre="$dir_MAGinator/maginator/workflow/prescreening_genes.Snakefile"
# snake_sign="$dir_MAGinator/maginator/workflow/signature_genes.Snakefile"

# dir_data="/home/projects/dtu_00009/people/henspi/git/Screener/data/simulated_data/MAGinator"

# echo "$snake_pre"
# # part of config reads=reads.csv contigs=contigs.fasta vamb=clusters.tsv params=test/parameters.tab
# snakemake --use-conda -s "$snake_pre" --resources mem_gb=180 --config wd="$dir_data"\
#  contigs="$dir_data/unused/contigs.fasta" vamb="$dir_data/unused/clusters.tsv" params="$dir_data/unused/parameters.tab"\
#   --cores 39 --printshellcmds format_conversion prescreening_genes
# #reads=reads.csv 
# echo "$snake_sign"
# snakemake --use-conda -s "$snake_sign" --resources mem_gb=180 \
# --config wd="$dir_data" n_refined_genes=500 contigs="$dir_data/unused/contigs.fasta" vamb="$dir_data/unused/clusters.tsv" params="$dir_data/unused/parameters.tab"\
#  --cores 39 --printshellcmds #--until gene_counts

