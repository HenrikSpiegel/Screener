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

import numpy as np
import configparser
from pathlib import Path

### Frontmatter
CONFIG_FILE = Path('config/analysis_ibdmdb.ini')
assert(Path(CONFIG_FILE).is_file())

config = configparser.ConfigParser()
config.read(CONFIG_FILE)

WD_DATA = Path(config.get("ProjectWide","WorkingDirData"))
WD_DATA.mkdir(exist_ok=True)
print(WD_DATA)
LOGLEVEL = config.get("ProjectWide","LoggingLevel")

#### Prepare datasources: This is somewhat specific for the analysis.
antismash_dirs = config.get("DataSources", "AntismashDirs")

sample_read_dir = Path(config.get("DataSources", "ReadPreprocessed"))

## Create symlink to get expected structure.
clean_samples = list(sample_read_dir.glob("clean_*"))
sample_names = list(set(fp.name.split("_")[1] for fp in clean_samples))

dir_preprocessed = WD_DATA/"preprocessed"
dir_preprocessed.mkdir(exist_ok=True)
for sample_name in list(sample_names):
    sym_dir = dir_preprocessed/sample_name
    sym_dir.mkdir(exist_ok=True)
    for file_name in [f"clean_{sample_name}_{i}.fastq" for i in [1,2]]:
        if (sym_dir/file_name).is_symlink():
            continue
        (sym_dir/file_name).symlink_to(sample_read_dir/file_name)



# Add tasks and dependencies:
job_id_map = dict()
dependencies = []


## Load the antismash called regions into a single file - is included as a step as it can be memory intensive.
dependencies.append(
    ("antismash_collect", )
)
job_id_map["antismash_collect"] = AddToQue(
    command=f"""\
python scripts/extract_antismash_ibdmdb.py\
 --antismash-dirs {antismash_dirs}\
 --outdir {WD_DATA}/antismash\
 --threads 10\
""",
    name="antismash_extract",
    success_file=WD_DATA / "antismash/.success"
)

# dependencies.append(
#     ("preprocess_describe",)
#     )
# job_id_map["preprocess_describe"] = AddToQue(
#     command=f"""\
# python scripts/extract_average_readlength.py\
#  --top-dir {WD_DATA}/preprocessed\
#  --fuzzy-dataset-names 'SR*'\
#  -o {WD_DATA}/preprocessed/average_lenghts.json\
#  --threads 8\
# """,
#     name="preprocess.describe",
#     success_file=WD_DATA / "preprocessed/.success_describe",
#     loglvl=LOGLEVEL
# )


## add pairwise blast of found bgcs.
dependencies.append(
    ('antismash_collect', 'blast_pw')
)
job_id_map['blast_pw'] = PairwiseBlast(
    fasta_file = WD_DATA / "antismash/combined_bgc.fa",
    output_dir=  WD_DATA / "blast_pairwise/input_bgc",
    loglvl=LOGLEVEL
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
    ("mcl_clustering", "analysis_12")
)
job_id_map["analysis_12"] = AddToQue(
    command=f"""\
python analysis/12_inputbgc_with_clustering.py\
 --blast {WD_DATA}/blast_pairwise/input_bgc/pairwise_table_symmetric.tsv\
 --mcl-cluster-dir {WD_DATA}/mcl_clustering\
 -o {WD_DATA}/results/12_inputbgc_with_clustering
""",
    success_file= WD_DATA / "results/12_inputbgc_with_clustering/.success"
)

# dependencies.append(
#     ({'antismash_collect','mcl_clustering'}, 'catalogue_generation')
# )
# job_id_map['catalogue_generation'] = AddToQue(
#     command=f"""\
# python -m lib.catalogue_assembler\
#  --bgcfasta {WD_DATA / 'antismash/combined_bgc.fa'}\
#  --families {WD_DATA / 'mcl_clustering/out.blast_result.mci.I40.json'}\
#  -o {WD_DATA / 'catalogues'}\
#  --max-catalogue-size 10000
# """,
#     success_file=WD_DATA / 'catalogues/success_catalogues',
#     name='catalogue_generation',
#     loglvl=LOGLEVEL
# )

kmerquant_labels = set()

#split jobs into subgroups to ease viewing. Purely cosmetic for log.graph
sample_dirs = list(dir_preprocessed.iterdir())
chunk_size = 10
sample_dirs_groups = [sample_dirs[x:x+chunk_size] for x in range(0, len(sample_dirs), chunk_size)]

for group_i, samle_dirs_i in enumerate(sample_dirs_groups):
    for sample_dir in samle_dirs_i:
        sample_label = sample_dir.name
        job_label = f"kmerQuant.{str(group_i).zfill(2)}.{sample_label}"
        kmerquant_labels.add(job_label)
        dependencies.append(({'antismash_collect','mcl_clustering'}, job_label))

        job_id_map[job_label] = QuantifierKmer(
            read_files = list(sample_dir.iterdir()),
            fp_catalogue = WD_DATA / "catalogues/catalogues",
            output_dir = WD_DATA / f"kmer_quantification/{sample_label}",
            kmer_size = config.getint("KmerQuantification","KmerLength"),
            loglvl=LOGLEVEL
        )


# add kmer count collection
dependencies.append((kmerquant_labels, "count.collect"))
count_fuzzy_path = WD_DATA / "kmer_quantification/SR*/counts/*.counted"
count_matrices_dir = WD_DATA / "kmer_quantification/count_matrices"
job_id_map["count.collect"] = AddToQue(
    command = f"python3 scripts/collect_count_data_analysis.py",
    success_file=count_matrices_dir/".success",
    name="count.collect",
    loglvl=LOGLEVEL
)

## Analyse distribution of counts.
dir_ana_13 = WD_DATA/ "results/13_count_distribution"
dir_ana_13.mkdir(parents=True, exist_ok=True)

dependencies.append(('count.collect', 'analysis_13'))
job_id_map["analysis_13"] = AddToQue(
    command = f"""\
python analysis/13_count_distribution.py\
 --dir-count-matrices {count_matrices_dir}\
 -o {dir_ana_13}/
""",
    name = "analysis_13",
    success_file=dir_ana_13/".success"
)

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
 --retries 50\
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
 --readlen-const 100\
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
 --antismash-dir {WD_DATA}/antismash\
 --outfile {fp_cluster_representatives}\
""",
    name="MapQuant.representatives",
    success_file=WD_DATA/"map_quantification/.success_clustrep"
)

# add map quantifications.
mapquant_labels = set()
for group_i, samle_dirs_i in enumerate(sample_dirs_groups):
    for sample_dir in samle_dirs_i:
        sample_name = sample_dir.name
        job_label = f"mapQuant.{group_i}.{sample_name}"
        mapquant_labels.add(job_label)

        dependencies.append(({"mcl_representatives"}, job_label))

        job_id_map[job_label] = QuantifierMap(
            reads = list(sample_dir.iterdir()),
            reference=fp_cluster_representatives,
            output_dir = WD_DATA / f"map_quantification/{sample_name}",
            loglvl=LOGLEVEL
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
