import configparser
import argparse
import os, sys
import logging
from pathlib import Path
import json

#Will run the analysis pipeline for the simulated data

from pipeline.general_functions import get_log, raise_w_log
from scripts.qsub_kmerquantifier import QuantifierKmer

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
    job_ids = args.dependencies
    log = get_log(log_name="RunAnalysis_"+jobtag, lvl=logging.DEBUG)

    config = configparser.ConfigParser()
    config.read("config/project_config.ini")

    dataset_dir = Path("data/simulated_data/preprocessed") / jobtag
    sample_set = list(dataset_dir.glob("sample_*"))
    log.info(f"Found ({len(sample_set)}) samples in dataset")


    if not args.runQuantifierMap:
        log.info("Skipping QuantifierMap, missing flag -> --runQuantifierMap")
    else:
        pass


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
            sample_reads = list(sample_dir.glob("*reads.fq.gz"))
            log.debug(f"Found read files: {sample_reads}")

            api_kmer = QuantifierKmer(read_files=sample_reads, fp_catalogue=fp_catalogue, output_dir=sample_outdir, kmer_size = config.getint("KmerQuantification","KmerLength"))
            kmer_key = api_kmer.__class__.__name__ + sample_dir.name
            if api_kmer.successful():
                log.info(f"Already ran: {jobtag} {sample_dir.name} - skipping ...")
                continue
            api_kmer.preflight()
            api_kmer.set_qsub_args(jobtag=jobtag+sample_dir.name)
            api_kmer.add_to_que(test=False)
            job_ids[api_kmer] = api_kmer.job_id
            
