
import os, sys, time
import re
from typing import Callable
import pandas as pd
import io
import pathlib
import glob

from scripts.functions import submit2
def gen_prefix(reads_gb:int) -> str:
    return str(reads_gb).replace(".","_")+"GB"

def gen_qsub_args(working_dir:str=None, job_tag:str="", **kwargs):
    if working_dir is None:
        working_dir=os.getcwd()

    scriptname = os.path.basename(__file__).split(".")[0]
    scriptname = "_".join([scriptname, job_tag])
    qsub_args = dict(
        directory = working_dir,
        modules = "tools anaconda3/2021.11 perl samtools/1.13 camisim/1.3",
        runtime = 120,
        cores = 30,
        ram=100,
        group="dtu_00009",
        jobname=scriptname,
        output = os.path.join(working_dir, "logs", scriptname+ "_stdout"),
        error = os.path.join(working_dir, "logs", scriptname+ "_stderr")
    )
    qsub_args.update(kwargs)
    return qsub_args


def generate_syscall(reads_GB:int, working_dir:str =None) -> str:
    if working_dir is None:
        working_dir = os.getcwd()
    file_prefix = gen_prefix(reads_GB)
    config_fp = f"data/simulated_data/camisim/{file_prefix}_config.ini"
    call_run_camisim = f"""\
python /services/tools/camisim/1.3/metagenomesimulation.py {os.path.join(working_dir, config_fp)}
mv {os.path.join(working_dir, f"data/simulated_data/camisim/{file_prefix}/*sample_0")} {os.path.join(working_dir, f"data/simulated_data/camisim/{file_prefix}/sample_0")}
"""
    return call_run_camisim

def _generate_supportingfiles(reads_GB: float):
    """
    Generates the config files needed to run CAMISIM for frozen set of input genomes with variable mass of
    reads generated.

    """
    file_prefix = gen_prefix(reads_GB)
    working_dir = os.getcwd()
    id_to_genome_fh = f"data/simulated_data/camisim/id_map.tsv"
    if not os.path.isfile(id_to_genome_fh):
        id_to_genome_str = f"""\
Genome1,{working_dir}/data/simulated_data/input_genomes/NC_014328_1.fa
Genome2,{working_dir}/data/simulated_data/input_genomes/NZ_CP020566_1.fa
Genome3,{working_dir}/data/simulated_data/input_genomes/NZ_CP053893_1.fa
Genome4,{working_dir}/data/simulated_data/input_genomes/NZ_LT906445_1.fa
Genome5,{working_dir}/data/simulated_data/input_genomes/NZ_LT906470_1.fa
"""
        pd.read_csv(io.StringIO(id_to_genome_str), sep=",").to_csv(id_to_genome_fh, sep="\t", index=False)
   
    metadata_fh = f"data/simulated_data/camisim/meta.tsv"
    if not os.path.isfile(metadata_fh):
        metadata_str = """\
genome_ID,OTU,NCBI_ID,novelty_category
Genome1,x,748727,known_species
Genome2,x,39777,known_species
Genome3,x,1520,known_species
Genome4,x,29466,known_species
Genome5,x,248315,known_species
    """
        pd.read_csv(io.StringIO(metadata_str), sep=",").to_csv(metadata_fh, sep="\t", index=False)
    distribution_fh = "data/simulated_data/camisim/distribution.tsv"
    if not os.path.isfile(distribution_fh):
        distribution_str = """\
Genome1,0.2
Genome2,0.2
Genome3,0.2
Genome4,0.2
Genome5,0.2
"""
        pd.read_csv(io.StringIO(distribution_str), sep=",").to_csv(distribution_fh, sep="\t", index=False)

def _generate_config(reads_GB: float, seed=20220709):
    file_prefix = gen_prefix(reads_GB)
    working_dir = os.getcwd()
    outdir = os.path.join(working_dir, f"data/simulated_data/camisim/{file_prefix}")

    #id_to_genome_fh = os.path.join(working_dir, f"data/simulated_data/camisim/{file_prefix}_id_map.tsv")
    #metadata_fh = os.path.join(working_dir, f"data/simulated_data/camisim/{file_prefix}_meta.tsv")
    config_fp = f"data/simulated_data/camisim/{file_prefix}_config.ini"
    config_content = f"""\
[Main]
# maximum number of processes
max_processors=30

# 0: community design + read simulator,
# 1: read simulator only
phase=0

# ouput directory, where the output will be stored (will be overwritten if set in from_profile)
output_directory={outdir}

# temporary directory
temp_directory=/tmp

# gold standard assembly
gsa=False

# gold standard for all samples combined
pooled_gsa=False

# anonymize sequences?
anonymous=True

# compress data (levels 0-9, recommended is 1 the gain of higher levels is not too high)
compress=1

# id of dataset, used in foldernames and is prefix in anonymous sequences
dataset_id={file_prefix}

# Read Simulation settings, relevant also for from_profile
[ReadSimulator]
# which readsimulator to use:
#           Choice of 'art', 'wgsim', 'nanosim', 'pbsim'
type=art

# Samtools (http://www.htslib.org/) takes care of sam/bam files. Version 1.0 or higher required!
# file path to executable
samtools=/services/tools/camisim/1.3/tools/samtools-1.3/samtools

# file path to read simulation executable
readsim=/services/tools/camisim/1.3/tools/art_illumina-2.3.6/art_illumina

#error profiles:
#for ART:
#HiSeq 150bp: hi150
#MBARC-26 150bp: mbarc
#custom profile (see below): own
#for wgsim:
#error rate as <float> (e.g. 0.05 for 5% error rate)
#blank for nanosim and wgsim
profile=mbarc

# Directory containing error profiles (can be blank for wgsim)
error_profiles=/services/tools/camisim/1.3/tools/art_illumina-2.3.6/profiles/

#paired end read, insert size (not applicable for nanosim)
fragments_size_mean=270
fragment_size_standard_deviation=27

# Only relevant if not from_profile is run:
[CommunityDesign]
# specify the samples size in Giga base pairs
size={reads_GB}

distribution_file_paths= data/simulated_data/camisim/distribution.tsv

# how many different samples?
number_of_samples=1

# how many communities
num_communities=1

# directory containing the taxdump of ncbi, version from 22.02.2017 is shipped
# "nodes.dmp"
# "merged.dmp"
# "names.dmp"
ncbi_taxdump=/services/tools/camisim/1.3/tools/ncbi-taxonomy_20170222.tar.gz

# the strain simulator for de novo strain creation
strain_simulation_template=/services/tools/camisim/1.3/scripts/StrainSimulationWrapper/sgEvolver/simulation_dir/

# define communities: [community<integer>]
[community0]
# information about all included genomes:
# can be used for multiple samples
metadata=data/simulated_data/camisim/meta.tsv
id_to_genome_file=data/simulated_data/camisim/id_map.tsv

# how many genomes do you want to sample over all?
genomes_total=5
num_real_genomes=5

# how many genomes per species taxon
#   (species taxon will be replaced by OTU-cluster later on)
max_strains_per_otu=1
ratio=1

# which kind of different samples do you need?
#   replicates / timeseries_lognormal / timeseries_normal / differential
mode=

# Part: community design
# Set parameters of log-normal and normal distribution, number of samples
# sigma > 0; influences shape (higher sigma -> smaller peak and longer tail),
log_sigma=2

# mu (real number) is a parameter for the log-scale
log_mu=1

# do you want to see a distribution before you decide to use it? yes/no
view=no\
"""
    with open(config_fp, "w") as fh:
        fh.write(config_content)


def preflight(reads_GB, seed=20220709):
    ## should raise error if outdir is full
    ## outdir can be read from config file.
    file_prefix = gen_prefix(reads_GB)
    pathlib.Path("data/simulated_data/camisim").mkdir(parents=True, exist_ok=True)
    
    #check if the output dir already exists.
    outdir = f"data/simulated_data/camisim/{file_prefix}"
    if os.path.isdir(outdir):
        raise RuntimeError("Output directory already exists - camisim cancelled. -> "+outdir)

    _generate_supportingfiles(reads_GB=reads_GB)
    _generate_config(reads_GB=reads_GB)

def output_exists(reads_gb: int) -> bool:
    """Checks if the output from the run already exists.
    Can be used as a check to see if the job ran succesfully,
    or if the job doesn't need to run.
    """
    reads_fp =  f"data/simulated_data/camisim/{gen_prefix(reads_gb)}/sample_0/reads/anonymous_reads.fq.gz"
    return os.path.isfile(reads_fp)
    # reads_dir = f"data/simulated_data/camisim/{gen_prefix(reads_gb)}/sample_0/reads"
    # if not os.path.isdir(reads_dir):
    #     print("CAMISIM: Output dir not found")
    #     return False
    # file_names = os.listdir(reads_dir)
    # n_output_fasta_files = len(re.findall("(Genome[1-5][1,2])", "".join(file_names)))
    # print(f"CAMISIM: Found ({n_output_fasta_files}/10) outputfiles", file=sys.stderr)
    # return n_output_fasta_files == 10 #TODO: probably shouldn't be hardcoded.

if __name__ == "__main__":
    import argparse
    ## Front matter - handle input parsing
    parser = argparse.ArgumentParser()
    parser.add_argument('--readGB', required=True, help="Size out output(reads) in giga bases")
    args = parser.parse_args()


    jobtag   = gen_prefix(args.readGB)
    preflight(reads_GB=args.readGB)
    qsub_kwargs = gen_qsub_args(job_tag=jobtag)
  
    runid_camisim = submit2(command=generate_syscall(reads_GB=args.readGB), **qsub_kwargs, test=False)
    
    print("camisim jobid: " + runid_camisim, file=sys.stderr)



    



