import configparser
import argparse
import os, sys
import logging
from pathlib import Path
import json

#Will run the analysis pipeline for the simulated data

from pipeline.general_functions import get_log, raise_w_log
from scripts.qsub_kmerquantifier import QuantifierKmer
from scripts.qsub_mapquantifier import QuantifierMap

if __name__ == "__main__":
    desc = """\
Pipeline for generating and preprocessing n samples for a single value of readsGB.
The number of samples are defined in config/project_config.ini.
    """
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("--readsGB", required=True, help="VolumeOfReads (float)")
    parser.add_argument("--dependencies", default='{}', type=json.loads, help="Dictionary of upstream qsub dependencies.(json.dumps format)")
    parser.add_argument("--runQuantifierMap", action="store_true", help="If flagged runs quantification-by-mapping pipeline")
    parser.add_argument("--runQuantifierKmer", action="store_true", help="If flagged runs quantification-by-kmer pipeline")
    args = parser.parse_args()

    jobtag = QuantifierKmer.gen_prefix(args.readsGB)
    job_ids = {}
    log = get_log(log_name="RunQuantification_"+jobtag, lvl=logging.INFO)

    config = configparser.ConfigParser()
    config.read("config/project_config.ini")

    dataset_dir = Path("data/simulated_data/camisim") / jobtag
 
    if args.dependencies:
        upstream_dependencies=list(args.dependencies.values())
        log.info(f"Found upstream dependencies: {upstream_dependencies}")
        log.info(f"Generating jobs on expected outputs matching sample count from config.")
        sample_set = [dataset_dir / f"sample_{x}" for x in range(config.getint("Simulation", "SimulatedSamples"))]
        log.info(f"Expecting ({len(sample_set)}) samples for data set.")
    else:
        upstream_dependencies=[]
        sample_set = list(dataset_dir.glob("sample_*"))
        log.info(f"Found ({len(sample_set)}) samples in dataset")

    if not args.runQuantifierKmer:
        log.info("Skipping QuantifierKmer, missing flag -> --runQuantifierKmer")
    else:
        fp_catalogue = Path("data/simulated_data/quantification_kmer/bgc_catalogues") #TODO: use our own catalogue
        if not fp_catalogue.is_dir():
                log.warning("Generating catalogue files")
                os.system(f"python3 scripts/kmer_gen_catalogue.py --fastas data/simulated_data/antismash/input_genomes/combined_bgc.fa -o {fp_catalogue}")

        for sample_dir in sample_set:
            log.debug(f"Running for sample: {sample_dir.name}")
            sample_outdir = Path("data/simulated_data/quantification_kmer") / jobtag / sample_dir.name
            #sample_reads = list(sample_dir.glob("*reads.fq.gz"))
            sample_reads = list(sample_dir.glob("reads/anonymous_reads.fq.gz")) #circumventing preprocesss
            log.debug(f"Found read files: {sample_reads}")

            api_kmer = QuantifierKmer(read_files=sample_reads, fp_catalogue=fp_catalogue, output_dir=sample_outdir, kmer_size = config.getint("KmerQuantification","KmerLength"))
            kmer_key = api_kmer.__class__.__name__ + sample_dir.name
            if api_kmer.successful():
                api_kmer.log.info(f"Already ran: {jobtag} {sample_dir.name} - skipping ...")
                continue
            api_kmer.preflight()
            api_kmer.set_qsub_args(jobtag=jobtag+sample_dir.name, dependency=upstream_dependencies)
            api_kmer.add_to_que(test=False)
            job_ids[kmer_key] = api_kmer.job_id


    if not args.runQuantifierMap:
        log.info("Not running runQuantifierMap pipeline, set --runQuantifierFlag flag")
    else:
        reference = "data/simulated_data/antismash/input_genomes/combined_bgc.fa"
        outdir    = f"data/simulated_data/quantification_map/{jobtag}/"

        for sample_dir in sample_set:
            log.debug(f"Running for sample: {sample_dir.name}")
            sample_outdir = Path("data/simulated_data/quantification_map") / jobtag / sample_dir.name
            #sample_reads = list(sample_dir.glob("*reads.fq.gz"))
            sample_reads = list(sample_dir.glob("reads/anonymous_reads.fq.gz")) #circumventing preprocesss
            log.debug(f"Found read files: {sample_reads}")

            api_map = QuantifierMap(
                    reads       = sample_reads, 
                    reference   = reference, 
                    output_dir  = sample_outdir, 
                    minMapQ     = config.getint("MapQuantification", "MinMapQ"),
                    minBaseQ    = config.getint("MapQuantification", "MinBaseQ") ,
                    log         = log)
            if api_map.successful():
                api_map.log.info(f"Already ran: {jobtag} {sample_dir.name} - skipping ...")
                continue
            else:
                map_key = api_map.__class__.__name__ + sample_dir.name
                api_map.preflight(check_input=False)
                api_map.set_qsub_args(jobtag=jobtag+sample_dir.name, dependency=upstream_dependencies)
                api_map.generate_syscall() #not needed as we run the default
                api_map.add_to_que(test=False)
                job_ids[map_key] = api_map.job_id
    log.info(f"Added ({len(job_ids)}) quantification jobs")
    log.debug(f"Added quantification jobs: {json.dumps(job_ids)}")
    print(json.dumps(job_ids))
