import configparser
import argparse
import os, sys
import logging
import pathlib
import json


from pipeline.general_functions import get_log, raise_w_log, submit2
from scripts.qsub_camisim import Camisim
from scripts.qsub_preprocess import Preprocessor
from scripts.qsub_ncbiFetch import NCBIFetch
from scripts.qsub_antismash import Antismash

if __name__ == "__main__":
    desc = """\
Pipeline for generating and preprocessing n samples for a single value of readsGB.
The number of samples are defined in config/project_config.ini.
    """
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("--readsGB", required=True, help="VolumeOfReads (float)")
    parser.add_argument("--dependencies", default='{}', type=json.loads, help="Dictionary of upstream qsub dependencies.(json.dumps format)")
    args = parser.parse_args()

    jobtag = Camisim.gen_prefix(args.readsGB)
    job_ids = args.dependencies
    log = get_log(log_name="RunSimulation_"+jobtag, lvl=logging.DEBUG)

    config = configparser.ConfigParser()
    config.read("config/project_config.ini")

    ncbi_ids = config.get("Simulation", "GenomeIDs").strip("\n").split("\n")

    ### Fetching from ncbi ###

    fetch_outdir = "data/simulated_data/input_genomes"
    expected_output = [os.path.join(fetch_outdir,f"{x.replace('.', '_')}.fa") for x in ncbi_ids+["combined"]]
    api_fetch = NCBIFetch(ncbi_ids, fetch_outdir, log=log)
    fetch_key = api_fetch.__class__.__name__
    if fetch_key in job_ids:
        log.info(f"Fetch running in upstream dependency - skipping fetch ...")
        log.debug(f"k:{fetch_key}, v: {job_ids[fetch_key]}")
    elif all([os.path.isfile(x) for x in expected_output]): #TODO: move to api
        log.info(f"Input genomes found - skipping fetch ...")
    else: 
        api_fetch.preflight()
        api_fetch.set_qsub_args(jobtag=jobtag)
        api_fetch.add_to_que()
        job_ids[fetch_key] = api_fetch.job_id
    ###

    ### Antismash on input genomes ###

    input_genomes_combined = os.path.join(fetch_outdir, "combined.fa")
    antismash_outdir = "data/simulated_data/antismash/input_genomes"
    api_antismash = Antismash(fastafile=input_genomes_combined, outdir=antismash_outdir, log=log)
    antismash_key = api_antismash.__class__.__name__
    if antismash_key in job_ids:
        log.info(f"Antismash running in upstream dependency - skipping ...")
        log.debug(f"k:{antismash_key}, v: {job_ids[antismash_key]}")
    elif api_antismash.successful():
        log.info(f"Antismash on input genomes found - skipping ...")
    else: 
        api_antismash.preflight()
        api_antismash.set_qsub_args(jobtag=jobtag)
        api_antismash.add_to_que()
        job_ids[antismash_key] = api_antismash.job_id
    ###

    ### Camisim on input genomes ###
 
    camisim_main_dir = "data/simulated_data/camisim"
    
    api_camisim = Camisim(
        readsGB=args.readsGB,
        n_samples=config.getint("Simulation","SimulatedSamples"),
        outdir = camisim_main_dir,
        log=log
        )
    key_camisim = api_camisim.__class__.__name__
    if api_camisim.successful():
        log.info(f"Camisim output for {args.readsGB}GB found - skipping...")
    else:
        dependencies = [job_ids.get(x) for x in [fetch_key, antismash_key] if x in job_ids]
        api_camisim.preflight()
        api_camisim.set_qsub_args(jobtag=jobtag, dependency=dependencies)
        api_camisim.add_to_que()
        job_ids[key_camisim] = api_camisim.job_id
    ###

    ### Preprocess camisim ###
    preprocess_out = f"data/simulated_data/preprocessed/{jobtag}"

    read_files = [api_camisim.dir_datalabel / f"sample_{x}" / "reads" / "anonymous_reads.fq.gz" for x in range(config.getint("Simulation","SimulatedSamples"))]
    api_preprocessor = Preprocessor(reads_interleaved=read_files, outdir=preprocess_out)
    key_processor = api_preprocessor.__class__.__name__
    if api_preprocessor.successful():
        log.info(f"Preproccor output found - skipping ...")
    else:
        api_preprocessor.preflight()
        api_preprocessor.set_qsub_args(jobtag=jobtag, dependency=job_ids.get(key_camisim, []))
        api_preprocessor.add_to_que()
        job_ids[key_processor] = api_preprocessor.job_id
    ###
    
    out_ids = json.dumps(job_ids)
    log.info("Finished adding simulation jobs: "+out_ids)
    print(job_ids)