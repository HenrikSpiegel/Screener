import json
from qsub_modules.ncbiFetch import NCBIFetch
from qsub_modules.antismash import Antismash
from qsub_modules.camisim import Camisim
from qsub_modules.preprocess import Preprocessor
from qsub_modules.blastn_pw import PairwiseBlast
from qsub_modules.add_to_que import AddToQue
from qsub_modules.kmerquantifier import QuantifierKmer

from pipeline.pipeline_base import PipelineBase


import numpy as np
import configparser
from pathlib import Path
from Bio import SeqIO


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
camisim_labels = ['camisim.'+label for label in GB_LABELS]
dependencies.append(('fetch', set(camisim_labels)))
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
preprocess_labels   =  ['preprocess.'+label for label in GB_LABELS]
preprocess_output_directories  =  [f"data/simulated_data/preprocessed/{label}" for label in GB_LABELS]
dependencies.extend([(cam, pre) for cam, pre in zip(camisim_labels, preprocess_labels)])

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
     for label, input_files, outdir in zip(preprocess_labels, input_file_sets, preprocess_output_directories)
    }
)


## Analysis part

# Add demonstration of blast+similarity.



dependencies.append(
    ('antismash', 'blast_demo_prep')
)
job_id_map["blast_demo_prep"] = AddToQue(
    command="python scripts/blast_demo_prep.py --n_shuffled_sequences 2 --n_chunks 5",
    success_file="data/simulated_data/blast_pairwise/demo/.success_prep",
    name="blast_demo_prep"
)

blast_demo_dir = Path("data/simulated_data/blast_pairwise/demo")
blast_demo_file = blast_demo_dir/"demo.fa"
dependencies.append(
    ('blast_demo_prep', 'blast_demo')
)
job_id_map['blast_demo'] = PairwiseBlast(
    fasta_file = blast_demo_file,
    output_dir= blast_demo_dir
)
dependencies.append(
    ('blast_demo', 'analysis_07_demo')
)
job_id_map['analysis_07_demo'] = AddToQue(
    command=f"""\
python analysis/07_blast_visualisation.py\
 --fasta {blast_demo_file}\
 --blast data/simulated_data/blast_pairwise/demo/combined_blast_results.tsv\
 --similarity-table data/simulated_data/blast_pairwise/demo/pairwise_table_symmetric.tsv\
 -o results/07_blast_visualisation/demo\
 -ph 500\
 -pw 1000
""",
    success_file='results/07_blast_visualisation/demo/.success',
    name='analysis_07_demo',
)

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
    ('blast_pw', 'analysis_01')
)
job_id_map['analysis_01'] = AddToQue(
    command='python analysis/01_compare_input_bgc.py',
    success_file='results/01_compare_input_bgc/success',
    name='analysis_01',
)

# add blast visualization

dependencies.append(
    ('blast_pw', 'analysis_07')
)
job_id_map['analysis_07'] = AddToQue(
    command="""\
python analysis/07_blast_visualisation.py\
 --fasta data/simulated_data/antismash/input_genomes/combined_bgc.fa\
 --blast data/simulated_data/blast_pairwise/input_bgc/combined_blast_results.tsv\
 --similarity-table data/simulated_data/blast_pairwise/input_bgc/pairwise_table_symmetric.tsv\
 -ph 1800\
 -pw 1000
""",
    success_file='results/07_blast_visualisation/.success',
    name='analysis_07',
)

# add family generation based on (in future) blast_pw
dependencies.append(
    ('blast_pw', "bgc_family_gen")
)
job_id_map['bgc_family_gen'] = AddToQue(
    command='python scripts/generate_bgc_families.py -o data/simulated_data/catalogues/family_dump.json',
    success_file='data/simulated_data/catalogues/success_fam_gen',
    name='bgc_family_gen',
)

# add catalogue generation for each family.
dependencies.append(
    ({'antismash','bgc_family_gen'}, 'catalogue_generation')
)
job_id_map['catalogue_generation'] = AddToQue(
    command="""\
python -m lib.catalogue_assembler\
 --bgcfasta data/simulated_data/antismash/input_genomes/combined_bgc.fa\
 --families data/simulated_data/catalogues/family_dump.json\
 -o data/simulated_data/catalogues
""",
    success_file='data/simulated_data/catalogues/success_catalogues',
    name='catalogue_generation',
)

# add kmer quantification.
kmerquant_labels = set()
for label in GB_LABELS:
    for sample in range(N_SAMPLES):
        job_label = f"kmerQuant.{label}.{sample}"
        kmerquant_labels.add(job_label)
        datadir = Path(f"data/simulated_data/preprocessed/{label}/sample_{sample}")
        dependencies.append(({'preprocess.'+label, "catalogue_generation"}, job_label))

        job_id_map[job_label] = QuantifierKmer(
            read_files = [
                datadir/"trimmed.anonymous_reads.fq.gz",  
                datadir/"trimmed.singleanonymous_reads.fq.gz"],
            fp_catalogue = "data/simulated_data/catalogues/catalogues",
            output_dir = f"data/simulated_data/kmer_quantification/{label}/sample_{sample}",
            kmer_size = config.getint("KmerQuantification","KmerLength")
        )

# add kmer count collection
dependencies.append((kmerquant_labels, "count_collect"))
count_fuzzy_path = "data/simulated_data/kmer_quantification/*GB/sample_*/counts/*.counted"
count_matrices_dir = Path("data/simulated_data/kmer_quantification/count_matrices")
job_id_map["count_collect"] = AddToQue(
    command = f"python scripts/collect_count_matrices.py --fuzzy_path '{count_fuzzy_path}' -o {count_matrices_dir}",
    success_file=count_matrices_dir/".success",
    name="count_collect"
)


pipeline_simulate = PipelineBase(
    pipe_name="SimulateData",
    dependencies = dependencies,
    job_map = job_id_map,
    iteration_sleep=15,
    max_workers=10,
    rerun_downstream=True,
    testing=False
)

if __name__ == "__main__":

    pipeline_simulate.run_pipeline()

    pass

