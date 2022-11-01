from qsub_modules.ncbiFetch import NCBIFetch
from qsub_modules.antismash import Antismash
from qsub_modules.camisim import Camisim
from qsub_modules.preprocess import Preprocessor
from qsub_modules.blastn_pw import PairwiseBlast
from qsub_modules.add_to_que import AddToQue

from pipeline.pipeline_base import PipelineBase


import numpy as np
import configparser
from pathlib import Path


#### Front matter
config = configparser.ConfigParser()
config.read("config/project_config.ini")

if config.get("Simulation", "ReadsGBStep") != "":
    rstep = config.getfloat("Simulation", "ReadsGBStep")
    rmin = config.getfloat("Simulation", "ReadsGBMin")
    rmax = config.getfloat("Simulation", "ReadsGBMax")
    precision = config.get("Simulation", "ReadsGBStep").split(".")[1].__len__() #Høker significant digits.
    gbs_to_run = np.round(np.arange(rmin, rmax+rstep, rstep), precision).tolist()
else:
    GBS_TO_RUN = []

if config.get("Simulation", "ReadsGBExtra", fallback=None):
    gbs_extra = [float(x.strip()) for x in config.get("Simulation", "ReadsGBExtra").split()]
    GBS_TO_RUN.extend(gbs_extra)


def gen_prefix(reads_gb:float) -> str:
        # returns stringified float with GB suffix
        # 1 -> 1_0GB
        # 0.1 -> 0_1GB
        return str(float(reads_gb)).replace(".","_")+"GB"
GB_LABELS = [gen_prefix(gb) for gb in GBS_TO_RUN]

N_SAMPLES = config.getint("Simulation","SimulatedSamples")

#### Add tasks and dependencies:
job_id_map = dict()
dependencies = []


# add fetch
ncbi_ids = config.get("Simulation", "GenomeIDs").strip("\n").split("\n")
job_id_map["fetch"] = NCBIFetch(
    id_list      = ncbi_ids,
    outdir = "data/simulated_data/input_genomes"
)

# add antismash
dependencies.append(('fetch', 'antismash'))
job_id_map['antismash'] = Antismash(
    fastafile = "data/simulated_data/input_genomes/combined.fa",
    outdir    = "data/simulated_data/antismash/input_genomes"
)

# add camisim:
camisim_labels = ['camisim_'+label for label in GB_LABELS]
dependencies.append(('antismash', set(camisim_labels)))
job_id_map.update(
    {
    label: Camisim(
        readsGB=gb,
        n_samples=config.getint("Simulation","SimulatedSamples"),
        outdir = Path("data/simulated_data/camisim"))
    for label, gb in zip(camisim_labels, GBS_TO_RUN)
    }
)


# add preprocess:
preprocess_labels   =  ['preprocess_'+label for label in GB_LABELS]
dependencies.extend([(cam, pre) for cam, pre in zip(camisim_labels, preprocess_labels)])
output_directories  =  [f"data/simulated_data/preprocessed/{label}" for label in GB_LABELS]
input_file_sets = [
  [
      Path("data/simulated_data/camisim")/label/ f"sample_{sample_n}" / "reads" / "anonymous_reads.fq.gz" 
      for sample_n in range(N_SAMPLES)
  ] 
  for label in GB_LABELS
]      
job_id_map.update(
    {
    label: Preprocessor(
        reads_interleaved=input_files,
        outdir = outdir)
     for label, input_files, outdir in zip(preprocess_labels, input_file_sets, output_directories)
    }
)


## Analysis part
# add pairwise blast of bgcs.
dependencies.append(
    ('antismash', 'blast_pw')
)
job_id_map['blast_pw'] = PairwiseBlast(
    fasta_file = "data/simulated_data/antismash/input_genomes/combined_bgc.fa",
    output_dir= Path('data/simulated_data/blast_pairwise/input_bgc')
)

# add pw analysis of bgc
dependencies.append(
    ('blast_pw', '01_analysis')
)
job_id_map['01_analysis'] = AddToQue(
    command='python analysis/01_compare_input_bgc.py',
    success_file='results/01_compare_input_bgc/success',
    name='01_analysis',
)

# add pw analysis of partial copsag set.
dependencies.append(
    ('blast_pw_copsag',)
)
job_id_map["blast_pw_copsag"] = PairwiseBlast(
    fasta_file="/home/projects/dtu_00009/people/henspi/copsac_bgc.fa",
    output_dir="data/simulated_data/blast_pairwise/copsag"
)

pipeline_simulate = PipelineBase(
    pipe_name="SimulateData",
    dependencies = dependencies,
    job_map = job_id_map,
    iteration_sleep=30,
    testing=False
)

if __name__ == "__main__":
    pipeline_simulate.run_pipeline()
    pass
