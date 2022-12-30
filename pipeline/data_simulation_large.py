from pipeline.pipeline_base import PipelineBase

from qsub_modules.add_to_que import AddToQue
from qsub_modules.ncbiFetch import NCBIFetch_Assemblies
from qsub_modules.antismash import Antismash
from qsub_modules.camisim import Camisim
from qsub_modules.blastn_pw import PairwiseBlast
from qsub_modules.preprocess import Preprocessor
from qsub_modules.mcl_clustering import MCLClustering
from qsub_modules.kmerquantifier import QuantifierKmer
from qsub_modules.mapquantifier import QuantifierMap
from qsub_modules.maginator import MAGinator

#from qsub_modules.maginator import MAGinator

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

####################################################################
###################### Fetch / generate data #######################
####################################################################

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

## camisim summary df.
dependencies.append((set(camisim_labels), 'camisim.describe_runs'))
job_id_map["camisim.describe_runs"] = AddToQue(
    command=f"""\
python scripts/camisim_combine_descriptions.py\
 --camisim-overview-files {WD_DATA}/camisim/*GB/simulation_overview.csv\
 --camisim-config {WD_DATA /"camisim/configs"/ (GB_LABELS[0]+"_config.ini")} \
 --outfile {WD_DATA}/camisim/simulation_overview_full.tsv\
""",
    name="camisim.describe",
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

# summaries preprocess (get avg processed read-length)
dependencies.append(
    (set(preprocess_labels), "preprocess.describe")
    )
job_id_map["preprocess.describe"] = AddToQue(
    command=f"""\
python scripts/extract_average_readlength.py\
 --top-dir {WD_DATA}/preprocessed\
 --fuzzy-dataset-names '*GB/sample*'\
 -o {WD_DATA}/preprocessed/average_lenghts.json\
 --threads 8\
""",
    name="preprocess.describe",
    success_file=WD_DATA / "preprocessed/.success_describe",
    loglvl=LOGLEVEL
)


####################################################################
###################### Data analysis steps #########################
####################################################################

######## Cluster and generate catalogues

## add pairwise blast of found bgcs.
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
    ('blast_pw', 'analysis.01')
)
job_id_map['analysis.01'] = AddToQue(
    command=f"""\
python analysis/01_compare_input_bgc_genera.py\
 --blast {WD_DATA}/blast_pairwise/input_bgc/pairwise_table_symmetric.tsv\
 --genera-table-dir {WD_DATA}/input_genomes\
 -o {WD_DATA}/results/01_compare_input_bgc_genera\
""",        
    name='analysis.01',
    loglvl=LOGLEVEL,
    success_file= WD_DATA / "results/01_compare_input_bgc_genera/.success"
)

## Clustering of found bgcs.
dependencies.append(
    ('blast_pw', 'mcl_clustering')
)
job_id_map["mcl_clustering"] = MCLClustering(
    blast_file = WD_DATA / "blast_pairwise/input_bgc/combined_blast_results.tsv",
    output_dir = WD_DATA / "mcl_clustering",
    loglvl=LOGLEVEL,
)

# add analysis of clustering.
dependencies.append(
    ("mcl_clustering", "analysis.12")
)
job_id_map["analysis.12"] = AddToQue(
    command=f"""\
python analysis/12_inputbgc_with_clustering.py\
 --blast {WD_DATA}/blast_pairwise/input_bgc/pairwise_table_symmetric.tsv\
 --genera-table-dir {WD_DATA}/input_genomes\
 --mcl-cluster-dir {WD_DATA}/mcl_clustering\
 -o {WD_DATA}/results/12_inputbgc_with_clustering
""",
    success_file= WD_DATA / "results/12_inputbgc_with_clustering/.success"
)

### Generate catalogues from clusters.

dependencies.append(
    ({'antismash','mcl_clustering'}, 'catalogue_generation')
)
job_id_map['catalogue_generation'] = AddToQue(
    command=f"""\
python -m lib.catalogue_assembler\
 --bgcfasta {WD_DATA / 'antismash/input_genomes/combined_bgc.fa'}\
 --families {WD_DATA / 'mcl_clustering/out.blast_result.mci.I40.json'}\
 -o {WD_DATA / 'catalogues'}\
 --max-catalogue-size 10000
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
dependencies.append(({"count.collect", "preprocess.describe"}, "count.correct"))
job_id_map["count.correct"] = AddToQue(
    command = f"""\
python3 scripts/correct_count_matrices.py\
 --counts {count_matrices_dir}\
 --avg-readlengths {WD_DATA}/preprocessed/average_lenghts.json\
 -k {config.getint("KmerQuantification", "KmerLength")}\
 --est-read-err {config.getfloat("KmerQuantification", "PerBaseErrorRate")}\
 -o {count_corrected_dir}\
""",
    success_file=count_corrected_dir/".success",
    name="count.correct",
    loglvl=LOGLEVEL
)

## Analyse distribution of counts.
dir_ana_13 = WD_DATA/ "results/13_count_distribution"
dir_ana_13.mkdir(parents=True, exist_ok=True)

dependencies.append(('count.collect', 'analysis.13'))
job_id_map["analysis.13"] = AddToQue(
    command = f"""\
python analysis/13_count_distribution.py\
 --dir-count-matrices {count_matrices_dir}\
 -o {dir_ana_13}/
""",
    name = "analysis.13",
    success_file=dir_ana_13/".success"
)


# ### MAGINATOR stuff

# # prepare MAGinator input:
# dependencies.append(
#     ('count.collect', 'MAGinator.input_prep')
#     #({'count_collect', 'catalogue_generation'}, 'MAGinator.input_prep')
# )
# dir_MAGinator_top = WD_DATA / "MAGinator"
# job_id_map['MAGinator.input_prep'] = AddToQue(
#     command=f"""\
# python scripts/MAGinator_prepinput.py\
#  --catalogues {WD_DATA / 'catalogues/catalogues'}\
#  --count-matrix {count_matrices_dir/'counts_all.tsv'}\
#  --min-dataset 0.08\
#  -o {dir_MAGinator_top}\
# """,
#     success_file=dir_MAGinator_top/".success.input_prep",
#     name = 'MAGinator.input_prep',
#     loglvl=LOGLEVEL
# )

# # Run MAGinator (selected snakes)
# dependencies.append(
#     ('MAGinator.input_prep', 'MAGinator.main')
# )
# job_id_map['MAGinator.main'] = MAGinator(
#     MAGinator_dir="/home/projects/dtu_00009/people/henspi/git/MAGinator", #checkout of MAGinator repo.
#     MAGinator_wd=dir_MAGinator_top,
#     refined_set_size=config.getint("MAGinator", "RefinedSetSize"),
#     loglvl=LOGLEVEL
# )

# # Extract to flatfiles
# dependencies.append(
#     ('MAGinator.main', 'MAGinator.extract')
# )
# job_id_map['MAGinator.extract'] = AddToQue(
#     command=f"""
# mkdir -p {dir_MAGinator_top}/screened_flat
# Rscript --vanilla scripts/MAGinator_extract_results.R {dir_MAGinator_top}/collectionID_order.txt {dir_MAGinator_top}/signature_genes/screened {dir_MAGinator_top}/screened_flat
# """,
#     success_file=dir_MAGinator_top/".success_extract",
#     loglvl=LOGLEVEL
# )
# # We need to set the instance requirement in this specific way to away instance sharing between modules. 
# # We may be able to get around this by initing the qsub_requirements - but that is a lot of code to refactor.
# new_qsub_requirements = job_id_map['MAGinator.extract'].qsub_requirements.copy()
# new_qsub_requirements.update({'modules': 'tools gcc/7.4.0 intel/perflibs/2020_update4 R/4.0.0'})
# job_id_map['MAGinator.extract'].qsub_requirements = new_qsub_requirements

### MAGpy stuff
dependencies.append(
    ('count.collect', 'MAGpy.main')
)
job_id_map["MAGpy.main"] = AddToQue(
    command = f"""\
python -m lib.magpy\
 --count-files {WD_DATA}/kmer_quantification/count_matrices/counts_all.tsv\
 --meta-files {WD_DATA}/catalogues/metafiles/*.meta\
 --output {WD_DATA}/MAGpy\
 --full-output\
 --max-threads 20\
 --catalogue-size 500\
 --verbosity INFO\
 --rng-seed 2812\
 --step-sizes 40 35 30 25 20 15 10 5 2\
 --retries 15\
 --min-improvement 0.02\
""",
    name="MAGpy.main",
    success_file=WD_DATA/"MAGpy/.success"
)

dependencies.append(
    ('MAGpy.main', 'MAGpy.abundances')
)
job_id_map["MAGpy.abundances"] = AddToQue(
    command = f"""\
python scripts/counts_to_abundances.py\
 --count-matrices {count_matrices_dir}\
 --mag-screened {WD_DATA}/MAGpy/screened_flat\
 --read-len-json {WD_DATA}/preprocessed/average_lenghts.json\
 --kmer-len {config.get("KmerQuantification", "KmerLength")}\
 --error-rate-est {config.get("KmerQuantification", "PerBaseErrorRate")}\
 -o {WD_DATA}/abundances\
    """,
    name="to_abundances",
    success_file=WD_DATA/"abundances/.success"
)

############# 
## Quantification by mapping
fp_cluster_representatives = WD_DATA / "map_quantification/cluster_representatives.fa"
# Determine cluster refs:
dependencies.append(
    ('mcl_clustering', "mcl_representatives")
)
job_id_map["mcl_representatives"] = AddToQue(
    command=f"""\
python scripts/get_cluster_representatives.py\
 --family-json {WD_DATA}/mcl_clustering/out.blast_result.mci.I40.json\
 --similarity-matrix {WD_DATA}/blast_pairwise/input_bgc/pairwise_table_symmetric.tsv\
 --antismash-dir {WD_DATA}/antismash/input_genomes\
 --outfile {fp_cluster_representatives}\
""",
    name="MapQuant.representatives",
    success_file=WD_DATA/"map_quantification/.success_clustrep"
)

# add map quantifications.
mapquant_labels = set()
for label in GB_LABELS:
    for sample in range(N_SAMPLES):
        job_label = f"mapQuant.{label}.{sample}"
        mapquant_labels.add(job_label)
        datadir = WD_DATA / f"preprocessed/{label}/sample_{sample}"
        dependencies.append(({'preprocess.'+label, "mcl_representatives"}, job_label))

        job_id_map[job_label] = QuantifierMap(
            reads = [
                datadir / "trimmed.anonymous_reads.fq.gz",  
                datadir / "trimmed.singleanonymous_reads.fq.gz"],
            reference=fp_cluster_representatives,
            output_dir = WD_DATA / f"map_quantification/{label}/sample_{sample}",
            loglvl=LOGLEVEL
        )

######################################################################## 
######################## Result Investations ###########################

###### Analysis 5:
dir_ana_05 = WD_DATA / "results/05_simulation_results"
dir_ana_05.mkdir(parents=True, exist_ok=True)
dependencies.append(
    ("MAGpy.abundances", "analysis.05")
)
job_id_map["analysis.05"] = AddToQue(
    command=f"""\
python analysis/05_simulation_result_and_comparison.py\
 --simulation-overview {WD_DATA}/camisim/simulation_overview_full.tsv\
 --total-abundances {WD_DATA}/abundances/total_abundances.csv\
 --family-dump {WD_DATA}/mcl_clustering/out.blast_result.mci.I40.json\
 -o {dir_ana_05}\
""",
    name="analysis.05",
    success_file=dir_ana_05/".success"
)

###### Analysis 8
dir_ana_08 = WD_DATA / "results/08_mag_diagnostics"
dir_ana_08.mkdir(parents=True, exist_ok=True)
dependencies.append(
    ({"camisim.describe_runs",'MAGpy.main'}, 'analysis.08')
)
job_id_map['analysis.08'] = AddToQue(
    command=f"""\
python analysis/08_mag_diagnostics.py\
 --simulation-overview {WD_DATA}/camisim/simulation_overview_full.tsv\
 --family-json {WD_DATA}/mcl_clustering/out.blast_result.mci.I40.json\
 --count-mats {WD_DATA}/kmer_quantification/count_matrices_corrected\
 --mag-flat {WD_DATA}/MAGpy/screened_flat\
 -o {dir_ana_08}
""",
    name="analysis.08",
    success_file=dir_ana_08/".succes"
)

# ###### Analysis 9
# # #NOTE: Could potentially be split
# dir_ana_09 = WD_DATA / "results/09_mag_kmer_location"
# dir_ana_09_pileup = dir_ana_09/"pileup"
# dir_ana_09_pileup.mkdir(parents=True, exist_ok=True)
# dependencies.append(
#     ({"camisim.describe_runs", 'MAGinator.extract'}, 'analysis.09')
# )
# job_id_map['analysis.09'] = AddToQue(
#     command=f"""\
# #Run pileup
# module load samtools/1.14
# for ID_SAMPLE in $(cut -f1 {WD_DATA}/camisim/id_map.tsv)
# do
#     echo $ID_SAMPLE
#     samtools mpileup --min-MQ 0 --min-BQ 0 -a "{WD_DATA}/camisim/0_5GB/sample_0/bam/$ID_SAMPLE.bam"\
#      | awk '{{print $1","$2","$4}}' > "{dir_ana_09_pileup}/$ID_SAMPLE.csv"
# done
# module unload samtools/1.14

# python analysis/09_MAG_kmer_location.py\
#     --catalogues {WD_DATA}/catalogues/catalogues\
#     --family_dump {WD_DATA}/mcl_clustering/out.blast_result.mci.I40.json\
#     --mag-flat {WD_DATA}/MAGinator/screened_flat\
#     --antismash {WD_DATA}/antismash/input_genomes\
#     --simulation-overview {WD_DATA}/camisim/simulation_overview_full.tsv\
#     --count-matrices {WD_DATA}/kmer_quantification/count_matrices\
#     --camisim-id-map {WD_DATA}/camisim/id_map.tsv\
#     --pileup-dir {dir_ana_09_pileup}\
#     -o {dir_ana_09}\
# """,
#     success_file=dir_ana_09/".success",
#     loglvl=LOGLEVEL
# )

dependencies.append(
    (mapquant_labels.union({"MAGpy.abundances"}), "analysis.14")
)
job_id_map["analysis.14"] = AddToQue(
    command=f"""\
python analysis/14_compare_map_kmer.py\
 --dir-map-quant {WD_DATA}/map_quantification\
 --kmer-abundances {WD_DATA}/abundances/total_abundances.csv\
 --camisim-overview {WD_DATA}/camisim/simulation_overview_full.tsv\
 --cluster-json {WD_DATA}/mcl_clustering/out.blast_result.mci.I40.json\
 -o {WD_DATA}/results/14_compare_map_kmer\
""",
    success_file=WD_DATA/"results/14_compare_map_kmer.success",
    name="analysis.14"
)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--dry-run", action='store_true')
    parser.add_argument("-t", "--test-print", action='store_true')
    args = parser.parse_args()

    pipeline_simulate = PipelineBase(
        config_file=CONFIG_FILE,
        pipe_name=Path(__file__).stem,
        dependencies = dependencies,
        job_map = job_id_map,
        iteration_sleep=15,
        max_workers=12,
        rerun_downstream=True,
        testing=args.test_print
    )

    if args.dry_run:
        print("Dry-run nothing is added to the que.")
    else:
        pipeline_simulate.run_pipeline()


        #python -m qsub_modules.add_to_que --command 'python -m pipeline.data_simulation_large' --name pipeline --options runtime:4320 modules:"anaconda3/2021.11 graphviz/2.40.1" --success logs/rpipe.log.success