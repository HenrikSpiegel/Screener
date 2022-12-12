import shutil
from qsub_modules.base import Base

import logging
import os, sys, time
import re
from typing import Callable, List
import pandas as pd
import io
from pathlib import Path
import glob

class Camisim(Base):
    def __init__(self, fps_genome_overview:List[Path], fp_genomes_dir:Path, taxdump:Path=None, readsGB:float=0.2, n_samples:int=2, outdir:str="data/camisim", seed=29092022, log: logging.Logger=None, loglvl = "DEBUG") -> None:
        """
        Generates a synthetic metagenomic sample.
        The outdir determines the high level dir, while each dataset (readsGB) are
        placed in their own subdirectories.
        """
        self.fps_genome_overview = [Path(x) for x in fps_genome_overview]
        self.fp_genomes_dir = Path(fp_genomes_dir)

        if log:
            self.add_external_log(log)

        self.loglvl = loglvl
        self.readsGB = readsGB
        self.datalabel = self.gen_prefix(readsGB)
        self.n_samples = n_samples
        self.outdir = Path(outdir)
        self.seed = seed
        if not taxdump:
            taxdump="/services/tools/camisim/1.3/tools/ncbi-taxonomy_20170222.tar.gz"
            self.log.warning("Using CAMISIM packaged dump - may be outdated -> "+taxdump)
            
        self.taxdump = Path(taxdump)
        
        # We could include that the class runs using less hardcoded reference to input genomes.
        # Reference to genomes could be via config.ini.
        self.fp_meta    = self.outdir / "meta.tsv"   
        self.fp_id_map  = self.outdir / "id_map.tsv"
        self.fp_distri  = self.outdir / "distribution.tsv"

        self.dir_datalabel  = self.outdir / self.datalabel
        self.dir_temp = self.dir_datalabel / "tmp"
        
        self.fp_config      = self.outdir / 'configs' / f"{self.datalabel}_config.ini"

        self.success_file = self.dir_datalabel / ".success"
        

    qsub_requirements = dict(
        modules = "tools anaconda3/2021.11 perl samtools/1.13 camisim/1.3",
        runtime = 48*60, #With the larger genome sets we are pushing at least 6 hours at 0.5GB-5samples
        cores = 38,
        ram=188
        )

    def generate_supporting_files_from_df(self, df_camisim:pd.DataFrame):

        self.log.debug(f"writting -> {self.fp_id_map}")
        df_camisim[["genome_ID","fp"]].to_csv(self.fp_id_map, sep="\t", index=False, header=False)#, encoding='utf-8-sig')

        self.log.debug(f"writting -> {self.fp_meta}")
        df_camisim[["genome_ID","OTU","NCBI_ID","novelty_category"]].to_csv(self.fp_meta, sep="\t", index=False)#, encoding='utf-8-sig')

#     def generate_supporting_files(self):
#         """
#         Generates the config files needed to run CAMISIM for frozen set of input genomes with variable mass of
#         reads generated.
#         These could be generated running from the script allowing easier adjustments.

#         #Note this is done slightly convoluted due to errors writing \t files in a manual fashion.
#         """
#         if not self.fp_id_map.is_file():
#             self.log.debug(f"writting -> {self.fp_id_map}")
#             id_to_genome_str = f"""\
# Genome1,data/simulated_data/input_genomes/NC_014328_1.fa
# Genome2,data/simulated_data/input_genomes/NZ_CP020566_1.fa
# Genome3,data/simulated_data/input_genomes/NZ_CP053893_1.fa
# Genome4,data/simulated_data/input_genomes/NZ_LT906445_1.fa
# Genome5,data/simulated_data/input_genomes/NZ_LT906470_1.fa
# """
#             pd.read_csv(io.StringIO(id_to_genome_str), sep=",").to_csv(self.fp_id_map, sep="\t", index=False)
    
#         if not self.fp_meta.is_file():
#             self.log.debug(f"writting -> {self.fp_meta}")
#             metadata_str = """\
# genome_ID,OTU,NCBI_ID,novelty_category
# Genome1,x,748727,known_species
# Genome2,x,39777,known_species
# Genome3,x,1520,known_species
# Genome4,x,29466,known_species
# Genome5,x,248315,known_species
#     """
#             pd.read_csv(io.StringIO(metadata_str), sep=",").to_csv(self.fp_meta, sep="\t", index=False)

#         if not self.fp_distri.is_file():
#             self.log.debug(f"writting -> {self.fp_distri}")
#             distribution_str = """\
# Genome1,0.2
# Genome2,0.2
# Genome3,0.2
# Genome4,0.2
# Genome5,0.2
# """
#             pd.read_csv(io.StringIO(distribution_str), sep=",").to_csv(self.fp_distri, sep="\t", index=False)


    def generate_config(self):
        config_content = f"""\
[Main]
# maximum number of processes
max_processors={self.qsub_requirements["cores"]-1}
seed={self.seed}

# 0: community design + read simulator,
# 1: read simulator only
phase=0

# ouput directory, where the output will be stored (will be overwritten if set in from_profile)
output_directory={self.dir_datalabel}

# temporary directory #It appears that we run out of space on /tmp maybe we try on harddrive instead - will be slower.
temp_directory={self.dir_datalabel}/tmp 

# gold standard assembly
gsa=False

# gold standard for all samples combined
pooled_gsa=False

# anonymize sequences?
anonymous=True

# compress data (levels 0-9, recommended is 1 the gain of higher levels is not too high)
compress=1

# id of dataset, used in foldernames and is prefix in anonymous sequences
dataset_id={self.datalabel}

# Read Simulation settings, relevant also for from_profile
[ReadSimulator]
# which readsimulator to use:
#           Choice of 'art', 'wgsim', 'nanosim', 'pbsim'


samtools=/services/tools/camisim/1.3/tools/samtools-1.3/samtools

#Normal use
type=art
readsim=/services/tools/camisim/1.3/tools/art_illumina-2.3.6/art_illumina
profile=mbarc
error_profiles=/services/tools/camisim/1.3/tools/art_illumina-2.3.6/profiles/

#Testing specific error 
#type=wgsim
#readsim=/services/tools/camisim/1.3/tools/wgsim/wgsim
#profile=0.03
#error_profiles=

#paired end read, insert size (not applicable for nanosim)
fragments_size_mean=270
fragment_size_standard_deviation=27

# Only relevant if not from_profile is run:
[CommunityDesign]
# specify the samples size in Giga base pairs
size={self.readsGB}

#Note each sample must have its won distribution_file_path, even if they are the same.
#distribution_file_paths={self.fp_distri}
#distribution_file_paths={",".join(self.fp_distri.__str__() for x in range(self.n_samples))}

# how many different samples?
number_of_samples={self.n_samples}

# how many communities
num_communities=1

# directory containing the taxdump of ncbi, version from 22.02.2017 is shipped
# "nodes.dmp"
# "merged.dmp"
# "names.dmp"
#ncbi_taxdump=/services/tools/camisim/1.3/tools/ncbi-taxonomy_20170222.tar.gz
ncbi_taxdump={self.taxdump}

# the strain simulator for de novo strain creation
strain_simulation_template=/services/tools/camisim/1.3/scripts/StrainSimulationWrapper/sgEvolver/simulation_dir/

# define communities: [community<integer>]
[community0]
# information about all included genomes:
# can be used for multiple samples
metadata={self.fp_meta}
id_to_genome_file={self.fp_id_map}

# how many genomes do you want to sample over all?
genomes_total={self.n_genomes}
num_real_genomes={self.n_genomes}

# how many genomes per species taxon
#   (species taxon will be replaced by OTU-cluster later on)
max_strains_per_otu=1000
ratio=1

# which kind of different samples do you need?
#   replicates / timeseries_lognormal / timeseries_normal / differential
mode=differential

# Part: community design
# Set parameters of log-normal and normal distribution, number of samples
# sigma > 0; influences shape (higher sigma -> smaller peak and longer tail),
log_sigma=0.25

# mu (real number) is a parameter for the log-scale
log_mu=0

# do you want to see a distribution before you decide to use it? yes/no
view=no\
"""
        self.fp_config.write_text(config_content)

    def preflight(self, check_input=False) -> None:


        if self.dir_datalabel.is_dir():
            if self.success_file.is_file():
                raise RuntimeError("Camisim Already Ran successfully - Please reset directory.")
            self.log.warning(f"Removing files from previous failed runs -> {self.dir_datalabel}")
            shutil.rmtree(self.dir_datalabel)
        self.dir_temp.mkdir(parents=True, exist_ok=True)

        # if check_input:
        #     input_files = glob.glob("data/simulated_data/input_genomes/*.fa") #TODO: set to check from genomes from config.
        #     if len(input_files) != 6:
        #         msg = "Didnt find input files (5 genomes .fa + 1 combined.fa)"
        #         self.log.error(msg)
        #         raise IOError(msg)
        #     self.log.info("Found all input files")
        
        self.log.debug("Generating supporting files and run config.")
        config_folder = self.outdir / "configs"
        config_folder.mkdir(parents=True, exist_ok=True)

        # Grab specifications for genomes.
        df_assembly_overviews = pd.concat(
            pd.read_csv(fp, sep="\t").assign(genus = fp.stem.split("_")[0])
            for fp in self.fps_genome_overview
        )
        
        df_assembly_overviews.rename(columns= {'taxid':'NCBI_ID'}, inplace=True)
        given_columns = set(df_assembly_overviews.columns)
        required_columns = {"genus", "assembly_accession", "NCBI_ID"} #we could have chosen some softer names but this should work. NCBI_ID is taxid
        if (required_columns - given_columns) != set():
           raise RuntimeError(f"df_genome_specifications must contain: {required_columns}. Input missing -> {required_columns - given_columns}")
        df_camisim_pre =  df_assembly_overviews.loc[:,["genus", "assembly_accession", "NCBI_ID"]]
        df_camisim_pre["genome_ID"] = df_camisim_pre.genus + (df_camisim_pre.index+1).astype(str)
        df_camisim_pre.reset_index(drop=True, inplace=True)

        df_camisim_pre["OTU"] = "x"
        df_camisim_pre["novelty_category"] = "known_species"

        #df_camisim_pre["fp"] = [self.fp_genomes_dir/ (acc+".fna") for acc in df_camisim_pre.assembly_accession.values]
       
        #lets find the filepaths. We always try to unzip and rezip so we will remove .gz prefix if it exists.
        filepaths = []
        for assembly in df_camisim_pre.assembly_accession:
            assembly_fp = list(self.fp_genomes_dir.glob(assembly+"*"))
            if not len(assembly_fp) == 1:
                raise RuntimeError(f"Did not find 1 unique fp found ({len(assembly_fp)}) -> {assembly}")
            filepaths.append(assembly_fp[0].as_posix().replace(".gz",""))
        df_camisim_pre["fp"] = filepaths

        self.n_genomes = len(df_camisim_pre)
        self.log.info(f"Found ({self.n_genomes}) genomes to simulate")

        #Ensure taxdump exists.
        if not self.taxdump.exists():
            raise FileExistsError(f"Couldnt locate taxdump: -> {self.taxdump}")
        
        self.generate_supporting_files_from_df(df_camisim_pre)
        self.generate_config()

    def generate_syscall(self) -> None:
        # sets mem / cpu based on default qsub args.

        syscall=f"""\
#unzip files.
#echo "unzipping {self.fp_genomes_dir}"
#gunzip {self.fp_genomes_dir}/*.gz

python /services/tools/camisim/1.3/metagenomesimulation.py {self.fp_config}

#echo "zipping {self.fp_genomes_dir}"
#gzip {self.fp_genomes_dir}/*

#Cleanup generated sample names.
START=0
END={self.n_samples-1}
for (( c=$START; c<=$END; c++ ))
do
    cmd="mv {self.dir_datalabel}/*sample_$c {self.dir_datalabel}/sample_$c"
    $cmd
done

python scripts/camisim_describe_run.py -d {self.dir_datalabel}

"""    
        self._syscall = syscall
        return

    # def successful(self):
    #     success_file = self.dir_datalabel / "success"
    #     return success_file.exists()

    # @staticmethod
    # def is_success(dir_datalabel) -> str:
    #     return os.path.isfile(os.path.join(dir_datalabel, "success"))

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--readsGB', required=True, help="Size out output(reads) in giga bases")
    parser.add_argument('--nSamples', default=2, type=int, help="Number of sample to generate [2]")
    parser.add_argument('-o', default="data/simulated_data/camisim", help="Top directory, subdirectory for readGB is created[data/simulated_data/camisim]")
    args = parser.parse_args()

    api = Camisim(readsGB=args.readsGB, n_samples=args.nSamples, outdir=args.o)
    api.preflight(check_input=True)
    api.set_qsub_args(jobtag="test")
    api.generate_syscall() #not needed as we run the default
    api.add_to_que(test=False)