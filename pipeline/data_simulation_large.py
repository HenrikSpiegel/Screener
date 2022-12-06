from pipeline.pipeline_base import PipelineBase

from qsub_modules.add_to_que import AddToQue
from qsub_modules.ncbiFetch import NCBIFetch_Assemblies
from qsub_modules.antismash import Antismash
from qsub_modules.camisim import Camisim
from qsub_modules.blastn_pw import PairwiseBlast
from qsub_modules.preprocess import Preprocessor

import numpy as np
import configparser
from pathlib import Path
from Bio import SeqIO
import sys

### Frontmatter

CONFIG_FILE = Path('config/project_config_large.ini')
assert(Path(CONFIG_FILE).is_file())

config = configparser.ConfigParser()
config.read(CONFIG_FILE)

WD_DATA = Path(config.get("ProjectWide","WorkingDirData"))
print(WD_DATA)
LOGLEVEL = config.get("ProjectWide","LoggingLevel")

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

GENERA = config.get("Simulation", "Genera").strip("\n").split("\n")
GENERA_CLEANED = [x.replace(" ", "_") for x in GENERA]

#################### MAIN MATTER BUILD PIPELINE ####################
# Add tasks and dependencies:
job_id_map = dict()
dependencies = []

#### Fetch input genomes.
genera_selected = config.get("Simulation", "Genera").strip("\n").split("\n")
genera_table_fps = [WD_DATA / "input_genomes/genera_tables" / fn for fn in genera_selected]
print( genera_table_fps)

dependencies.append(
    ('fetch',)
)
job_id_map["fetch"] = NCBIFetch_Assemblies(
    genera_names     = genera_selected,
    genera_table_dir = WD_DATA / "input_genomes/genera_tables",
    output_dir       = WD_DATA / "input_genomes",
    loglvl           = LOGLEVEL
)

#### Antismash on input genomes.
dependencies.append(
    ('fetch', 'antismash')
)
job_id_map["antismash"] = Antismash(
    fastafile=WD_DATA / "input_genomes/combined_genomes.fna.gz",
    outdir= WD_DATA / "antismash/input_genomes"
)
# The job is much larger that initial expected - so lets give it more power
job_id_map["antismash"].qsub_requirements.update(
    dict(
        runtime = 24*60,
        cores = 38,
        ram = 180,
    )
)

### CAMISIM Simulations.
camisim_labels = ['camisim.'+label for label in GB_LABELS]
genome_overview_fps = [WD_DATA / "input_genomes" / (g+"_selected.tsv") for g in GENERA_CLEANED]
dependencies.append(('fetch', set(camisim_labels)))
job_id_map.update(
    {
    label: Camisim(
        fps_genome_overview = genome_overview_fps,
        fp_genomes_dir = WD_DATA / "input_genomes/genomes",
        taxdump=WD_DATA / "input_genomes/taxdump/ncbi_tax_dump.tar.gz",
        readsGB=gb,
        n_samples=config.getint("Simulation","SimulatedSamples"),
        outdir = WD_DATA / "camisim",
        loglvl=LOGLEVEL)
    for label, gb in zip(camisim_labels, GBS_TO_RUN)
    }
)

# add preprocess:
preprocess_labels   =  ['preprocess.'+label for label in GB_LABELS]
preprocess_output_directories  =  [WD_DATA / f"preprocessed/{label}" for label in GB_LABELS]
dependencies.extend([(cam, pre) for cam, pre in zip(camisim_labels, preprocess_labels)])

input_file_sets = [
  [
      WD_DATA / "camisim" /label/ f"sample_{sample_n}" / "reads" / "anonymous_reads.fq.gz" 
      for sample_n in range(N_SAMPLES)
  ] 
  for label in GB_LABELS
]      
job_id_map.update(
    {
    label: Preprocessor(
        reads_interleaved=input_files,
        outdir = outdir,
        loglvl=LOGLEVEL)
     for label, input_files, outdir in zip(preprocess_labels, input_file_sets, preprocess_output_directories)
    }
)


# # add pairwise blast of found bgcs.
# dependencies.append(
#     ('antismash', 'blast_pw')
# )
# job_id_map['blast_pw'] = PairwiseBlast(
#     fasta_file = WD_DATA / "antismash/input_genomes/combined_bgc.fa",
#     output_dir=  WD_DATA / "blast_pairwise/input_bgc",
#     loglvl=LOGLEVEL
# )

# # add pw analysis of bgc
# dependencies.append(
#     ('blast_pw', 'analysis_01')
# )
# job_id_map['analysis_01'] = AddToQue(
#     command=f"python analysis/01_compare_input_bgc.py --antismashdir {WD_DATA / 'antismash/input_genomes'} --blast-table {WD_DATA / 'blast_pairwise/input_bgc/pairwise_table_symmetric.tsv'} -o {WD_DATA /'results/01_compare_input_bgc'}",
#     success_file=WD_DATA /'results/01_compare_input_bgc/.success',
#     name='analysis_01',
#     loglvl=LOGLEVEL
# )

pipeline_simulate = PipelineBase(
    config_file=CONFIG_FILE,
    pipe_name=Path(__file__).stem,
    dependencies = dependencies,
    job_map = job_id_map,
    iteration_sleep=15,
    max_workers=12,
    rerun_downstream=True,
    testing=False
)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--dry-run", action='store_true')
    args = parser.parse_args()
    if args.dry_run:
        print("Dry-run nothing is added to the que.")
        pass
    else:
        pipeline_simulate.run_pipeline()