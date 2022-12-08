import json
from qsub_modules.ncbiFetch import NCBIFetch
from qsub_modules.antismash import Antismash
from qsub_modules.camisim import Camisim
from qsub_modules.preprocess import Preprocessor
from qsub_modules.blastn_pw import PairwiseBlast
from qsub_modules.add_to_que import AddToQue
from qsub_modules.kmerquantifier import QuantifierKmer
from qsub_modules.mapquantifier import QuantifierMap
from qsub_modules.maginator import MAGinator

from pipeline.pipeline_base import PipelineBase


import numpy as np
import configparser
from pathlib import Path
from Bio import SeqIO


#### Front matter

CONFIG_FILE = Path('config/project_config_init.ini')

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

#### Add tasks and dependencies:
job_id_map = dict()
dependencies = []


# add fetch
ncbi_ids, tax_ids = zip(*[x.split(",") for x in config.get("Simulation", "GenomeIDs").strip("\n").split("\n")])


job_id_map["fetch"] = NCBIFetch(
    acc_ids=ncbi_ids,
    tax_ids=tax_ids,
    outdir = WD_DATA / "input_genomes",
    loglvl=LOGLEVEL
)

# add antismash
dependencies.append(('fetch', 'antismash'))
job_id_map['antismash'] = Antismash(
    fastafile = WD_DATA / "input_genomes/combined.fa",
    outdir    = WD_DATA / "antismash/input_genomes",
    loglvl    = LOGLEVEL
)

# add camisim:
camisim_labels = ['camisim.'+label for label in GB_LABELS]
dependencies.append(('fetch', set(camisim_labels)))
job_id_map.update(
    {
    label: Camisim(
        fps_genome_overview=[WD_DATA / "input_genomes/ncbi_selected.tsv"],
        fp_genomes_dir=WD_DATA / "input_genomes",
        readsGB=gb,
        n_samples=config.getint("Simulation","SimulatedSamples"),
        outdir = WD_DATA / "camisim",
        loglvl=LOGLEVEL)
    for label, gb in zip(camisim_labels, GBS_TO_RUN)
    }
)

## camisim summary df.
dependencies.append((set(camisim_labels), 'camisim.collect'))
job_id_map["camisim.collect"] = AddToQue(
    command=f"""\
python scripts/camisim_combine_descriptions.py\
 --camisim-overview-files {WD_DATA}/camisim/*GB/simulation_overview.csv\
 --camisim-config {WD_DATA /"camisim/configs"/ (GB_LABELS[0]+"_config.ini")} \
 --outfile {WD_DATA}/camisim/simulation_overview_full.tsv\
""",
    name="camisim.collect",
    success_file=WD_DATA / "camisim/.success_collect"
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


## Analysis part

# Add demonstration of blast+similarity.
blast_demo_dir = WD_DATA / "blast_pairwise/demo"

dependencies.append(
    ('antismash', 'blast_demo_prep')
)
job_id_map["blast_demo_prep"] = AddToQue(
    command=f"""\
python scripts/blast_demo_prep.py\
 --seq-file {WD_DATA / 'antismash/input_genomes/combined_bgc.fa'}\
 -o {blast_demo_dir}\
 --n_shuffled_sequences 2\
 --n_chunks 5\
""",
    success_file= WD_DATA / "blast_pairwise/demo/.success_prep",
    name="blast_demo_prep",
    loglvl=LOGLEVEL
)


blast_demo_file = blast_demo_dir/"demo.fa"
dependencies.append(
    ('blast_demo_prep', 'blast_demo')
)
job_id_map['blast_demo'] = PairwiseBlast(
    fasta_file = blast_demo_file,
    output_dir= blast_demo_dir,
    loglvl = LOGLEVEL
)

dependencies.append(
    ('blast_demo', 'analysis_07_demo')
)
job_id_map['analysis_07_demo'] = AddToQue(
    command=f"""\
python analysis/07_blast_visualisation.py\
 --fasta {blast_demo_file}\
 --blast {blast_demo_dir / 'combined_blast_results.tsv'}\
 --similarity-table {blast_demo_dir / 'pairwise_table_symmetric.tsv'}\
 -o {WD_DATA / 'results/07_blast_visualisation/demo'}\
 -ph 500\
 -pw 1000
""",
    success_file = WD_DATA / 'results/07_blast_visualisation/demo/.success',
    name='analysis_07_demo',
    loglvl=LOGLEVEL
)

# add pairwise blast of bgcs.
dependencies.append(
    ('antismash', 'blast_pw')
)
job_id_map['blast_pw'] = PairwiseBlast(
    fasta_file = WD_DATA / "antismash/input_genomes/combined_bgc.fa",
    output_dir=  WD_DATA / "blast_pairwise/input_bgc",
    loglvl=LOGLEVEL
)

# add pw analysis of bgc
dependencies.append(
    ('blast_pw', 'analysis_01')
)
job_id_map['analysis_01'] = AddToQue(
    command=f"""\
python analysis/01_compare_input_bgc.py\
 --antismashdir {WD_DATA / 'antismash/input_genomes'}\
 --blast-table {WD_DATA / 'blast_pairwise/input_bgc/pairwise_table_symmetric.tsv'}\
 -o {WD_DATA /'results/01_compare_input_bgc'}
 """,
    success_file=WD_DATA /'results/01_compare_input_bgc/.success',
    name='analysis_01',
    loglvl=LOGLEVEL

)

# add blast visualization

dependencies.append(
    ('blast_pw', 'analysis_07')
)
job_id_map['analysis_07'] = AddToQue(
    command=f"""\
python analysis/07_blast_visualisation.py\
 --fasta {WD_DATA / 'antismash/input_genomes/combined_bgc.fa'}\
 --blast {WD_DATA / 'blast_pairwise/input_bgc/combined_blast_results.tsv'}\
 --similarity-table {WD_DATA / 'blast_pairwise/input_bgc/pairwise_table_symmetric.tsv'}\
 -o {WD_DATA / "results/07_blast_visualisation"}\
 -ph 1800\
 -pw 1000
""",
    success_file=WD_DATA / 'results/07_blast_visualisation/.success',
    name='analysis_07',
    loglvl=LOGLEVEL
)

# add family generation based on (in future) blast_pw
dependencies.append(
    ('blast_pw', "bgc_family_gen")
)
job_id_map['bgc_family_gen'] = AddToQue(
    command=f"python scripts/generate_bgc_families.py -o {WD_DATA / 'catalogues/family_dump.json'}",
    success_file=WD_DATA / 'catalogues/success_fam_gen',
    name='bgc_family_gen',
    loglvl=LOGLEVEL
)

# add catalogue generation for each family.
dependencies.append(
    ({'antismash','bgc_family_gen'}, 'catalogue_generation')
)
job_id_map['catalogue_generation'] = AddToQue(
    command=f"""\
python -m lib.catalogue_assembler\
 --bgcfasta {WD_DATA / 'antismash/input_genomes/combined_bgc.fa'}\
 --families {WD_DATA / 'catalogues/family_dump.json'}\
 -o {WD_DATA / 'catalogues'}\
 --max-catalogue-size 5000
""",
    success_file=WD_DATA / 'catalogues/success_catalogues',
    name='catalogue_generation',
    loglvl=LOGLEVEL
)

# add kmer quantification.
kmerquant_labels = set()
for label in GB_LABELS:
    for sample in range(N_SAMPLES):
        job_label = f"kmerQuant.{label}.{sample}"
        kmerquant_labels.add(job_label)
        datadir = WD_DATA / f"preprocessed/{label}/sample_{sample}"
        dependencies.append(({'preprocess.'+label, "catalogue_generation"}, job_label))

        job_id_map[job_label] = QuantifierKmer(
            read_files = [
                datadir / "trimmed.anonymous_reads.fq.gz",  
                datadir / "trimmed.singleanonymous_reads.fq.gz"],
            fp_catalogue = WD_DATA / "catalogues/catalogues",
            output_dir = WD_DATA / f"kmer_quantification/{label}/sample_{sample}",
            kmer_size = config.getint("KmerQuantification","KmerLength"),
            loglvl=LOGLEVEL
        )

# add kmer count collection
dependencies.append((kmerquant_labels, "count.collect"))
count_fuzzy_path = WD_DATA / "kmer_quantification/*GB/sample_*/counts/*.counted"
count_matrices_dir = WD_DATA / "kmer_quantification/count_matrices"
job_id_map["count.collect"] = AddToQue(
    command = f"python3 scripts/collect_count_matrices.py --fuzzy_path '{count_fuzzy_path}' -o {count_matrices_dir}",
    success_file=count_matrices_dir/".success",
    name="count.collect",
    loglvl=LOGLEVEL
)

# add kmer count correction
count_corrected_dir = count_matrices_dir.with_name("count_matrices_corrected")
dependencies.append(("count.collect", "count.correct"))
job_id_map["count.correct"] = AddToQue(
    command = f"""\
python3 scripts/correct_count_matrices.py\
 --counts {count_matrices_dir}\
 --reads-dir {WD_DATA}/preprocessed\
 --fuzzy-dataset-names '*GB/sample_*'\
 -k {config.getint("KmerQuantification", "KmerLength")}\
 --est-read-err {config.getfloat("KmerQuantification", "PerBaseErrorRate")}\
 -o {count_corrected_dir}\
""",
    success_file=count_corrected_dir/".success",
    name="count.correct",
    loglvl=LOGLEVEL
)


# parser.add_argument("--counts", required=True, type=Path, help="dir for raw count matrices")
# parser.add_argument("--reads-dir", required=True, type=Path, help="Dir holding processed reads (or raw if not processing)")
# parser.add_argument("--fuzzy-dataset-names", type=str, default='*GB/sample_*', help="fuzzy names for datasets in reads_dir, assumes parent*/sample* style naming.['*GB/sample_*']")
# parser.add_argument("-k", default=21, type=int)
# parser.add_argument("--est-read-err", default=0.015, type=float, help="Estimated average per base read error.")

# parser.add_argument("-o", required=True, type=Path, help="outdir for corrected matrices.")

# # add mapping quantification
# map_quant_labels = set()

# mapquantifier_reference = WD_DATA / "/antismash/input_genomes/combined_bgc.fa"

# for label in GB_LABELS:
#     for sample in range(N_SAMPLES):
#         job_label = f"mapQuant.{label}.{sample}"
#         map_quant_labels.add(job_label)
#         datadir = Path(fWD_DATA / "/preprocessed/{label}/sample_{sample}")
#         dependencies.append(({'preprocess.'+label, "antismash"}, job_label))

#         job_id_map[job_label] = QuantifierMap(
#                     reads        = [
#                                     datadir/"trimmed.anonymous_reads.fq.gz",  
#                                     datadir/"trimmed.singleanonymous_reads.fq.gz"],
#                     reference   = mapquantifier_reference, 
#                     output_dir  = fWD_DATA / "/map_quantification/{label}/sample_{sample}", 
#                     minMapQ     = config.getint("MapQuantification", "MinMapQ"),
#                     minBaseQ    = config.getint("MapQuantification", "MinBaseQ")
#         )



###### MAGinator part ######

# prepare MAGinator input:
dependencies.append(
    ('count.correct', 'MAGinator.input_prep')
    #({'count_collect', 'catalogue_generation'}, 'MAGinator.input_prep')
)
dir_MAGinator_top = WD_DATA / "MAGinator"
job_id_map['MAGinator.input_prep'] = AddToQue(
    command=f"""\
python scripts/MAGinator_prepinput.py\
 --catalogues {WD_DATA / 'catalogues/catalogues'}\
 --count-matrix {count_corrected_dir/'counts_all.tsv'}\
 -o {dir_MAGinator_top}\
""",
    success_file=dir_MAGinator_top/".success.input_prep",
    name = 'MAGinator.input_prep',
    loglvl=LOGLEVEL
)

# Run MAGinator (selected snakes)
dependencies.append(
    ('MAGinator.input_prep', 'MAGinator.main')
)
job_id_map['MAGinator.main'] = MAGinator(
    MAGinator_dir="/home/projects/dtu_00009/people/henspi/git/MAGinator", #checkout of MAGinator repo.
    MAGinator_wd=dir_MAGinator_top,
    refined_set_size=config.getint("MAGinator", "RefinedSetSize"),
    loglvl=LOGLEVEL
)

# Extract to flatfiles
dependencies.append(
    ('MAGinator.main', 'MAGinator.extract')
)
job_id_map['MAGinator.extract'] = AddToQue(
    command=f"""
mkdir -p {dir_MAGinator_top}/screened_flat
Rscript --vanilla scripts/MAGinator_extract_results.R {dir_MAGinator_top}/collectionID_order.txt {dir_MAGinator_top}/signature_genes/screened {dir_MAGinator_top}/screened_flat
""",
    success_file=dir_MAGinator_top/".success_extract",
    loglvl=LOGLEVEL
)
# We need to set the instance requirement in this specific way to away instance sharing between modules. 
# We may be able to get around this by initing the qsub_requirements - but that is a lot of code to refactor.
new_qsub_requirements = job_id_map['MAGinator.extract'].qsub_requirements.copy()
new_qsub_requirements.update({'modules': 'tools gcc/7.4.0 intel/perflibs/2020_update4 R/4.0.0'})
job_id_map['MAGinator.extract'].qsub_requirements = new_qsub_requirements

### MAG analysis

dir_ana_08 = WD_DATA / "results/09_MAG_kmer_location"
dependencies.append(
    ({"camisim.collect",'MAGinator.extract'}, 'analysis.08')
)
job_id_map['analysis.08'] = AddToQue(
    command=f"""\
python analysis/08_mag_diagnostics.py\
 --simulation-overview {WD_DATA}/camisim/simulation_overview_full.tsv\
 --family-json {WD_DATA}/catalogues/family_dump.json\
 --count-mats {WD_DATA}/kmer_quantification/count_matrices\
 --mag-flat {WD_DATA}/MAGinator/screened_flat\
 -o {dir_ana_08}
""",
    name="analysis.08",
    success_file=dir_ana_08/".succes"
)
# --camisim-config {WD_DATA}/camisim/configs/0_5GB_config.ini\
# --fuzzy-simulation '{WD_DATA}/camisim/*GB/simulation_overview.csv'\

# Run analysis on MAGinator output.
dir_ana_09 = WD_DATA / "results/09_mag_kmer_location"
dir_ana_09_pileup = dir_ana_09/"pileup"
dir_ana_09_pileup.mkdir(parents=True, exist_ok=True)
dependencies.append(
    ({"camisim.collect", 'MAGinator.extract'}, 'analysis_09')
)
job_id_map['analysis_09'] = AddToQue(
    command=f"""\
#Run pileup
module load samtools/1.14
for ID_SAMPLE in $(cut -f1 {WD_DATA}/camisim/id_map.tsv)
do
    echo $ID_SAMPLE
    samtools mpileup --min-MQ 0 --min-BQ 0 -a "{WD_DATA}/camisim/0_5GB/sample_0/bam/$ID_SAMPLE.bam"\
     | awk '{{print $1","$2","$4}}' > "{dir_ana_09_pileup}/$ID_SAMPLE.csv"
done
module unload samtools/1.14

python analysis/09_MAG_kmer_location.py\
    --catalogues {WD_DATA}/catalogues/catalogues\
    --family_dump {WD_DATA}/catalogues/family_dump.json\
    --mag-flat {WD_DATA}/MAGinator/screened_flat\
    --antismash {WD_DATA}/antismash/input_genomes\
    --simulation-overview {WD_DATA}/camisim/simulation_overview_full.tsv\
    --count-matrices {WD_DATA}/kmer_quantification/count_matrices\
    --camisim-id-map {WD_DATA}/camisim/id_map.tsv\
    --pileup-dir {dir_ana_09_pileup}\
    -o {dir_ana_09}\
""",
    success_file=dir_ana_09/".success_extract",
    loglvl=LOGLEVEL
)



pipeline_simulate = PipelineBase(
    config_file=CONFIG_FILE,
    pipe_name=Path(__file__).stem,
    dependencies = dependencies,
    job_map = job_id_map,
    iteration_sleep=15,
    max_workers=12,
    rerun_downstream=True,
    testing=True
)

if __name__ == "__main__":
    #pass
    pipeline_simulate.run_pipeline()



