from qsub_modules.ncbiFetch import NCBIFetch
from qsub_modules.antismash import Antismash
from qsub_modules.camisim import Camisim
from qsub_modules.preprocess import Preprocessor

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

pipeline_simulate = PipelineBase(
    pipe_name="SimulateData",
    dependencies = dependencies,
    job_map = job_id_map,
    iteration_sleep=10,
    testing=False
)

if __name__ == "__main__":
    pipeline_simulate.start()

