import glob
import os, sys
import shutil
import re
import logging
import pathlib

from pipeline.general_functions import submit2

import scripts.qsub_camisim as camisim
import scripts.qsub_ncbiFetch as ncbifetch
import scripts.qsub_antismash as antismash
from scripts.qsub_base import Base
from scripts.qsub_assemble import Assembler
from scripts.qsub_mapquantifier import QuantifierMap
from scripts.qsub_kmerquantifier import QuantifierKmer

api_base = Base()

working_dir = api_base.working_dir
ncbi_id_file = "config/ids_simulation_genomes.txt"

def get_log(tag, logdir = "logs/pipeline", lvl=logging.DEBUG):
    log_name = os.path.basename(__file__).split(".")[0]+"_"+jobtag
    logger = logging.getLogger(log_name)
    
    #remove old handlers:
    while logger.handlers:
        logger.removeHandler(logger.handlers[0])
    
    logger.setLevel(lvl)
    F = "[%(asctime)s %(name)s:%(funcName)s]%(levelname)s: %(message)s"
    formatter = logging.Formatter(F, datefmt='%d-%b-%y %H:%M:%S')
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    stream_handler.setLevel(lvl)
    logger.addHandler(stream_handler)
    
    if logdir:
        fp_log = os.path.join(logdir, log_name)
        if os.path.isfile(fp_log): 
            os.unlink(fp_log)
            
        pathlib.Path(os.path.dirname(fp_log)).mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(fp_log)
        file_handler.setFormatter(formatter)
        file_handler.setLevel(lvl)
        logger.addHandler(file_handler)
        logger.debug("logfile at -> "+fp_log)
    return logger

def raise_w_log(log, exception:Exception, msg):
    log.error(msg)
    raise exception(msg)

if __name__ == "__main__":
    ### FRONT MATTER ###
    import argparse
    ## Front matter - handle input parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("--readsGB", required=True, help="VolumeOfReads (float)")
    parser.add_argument("--runAssembly", action="store_true", help="If flagged runs assembly of simulated reads")
    parser.add_argument("--runQuantifierMap", action="store_true", help="If flagged runs quantification-by-mapping pipeline")
    parser.add_argument("--runQuantifierKmer", action="store_true", help="If flagged runs quantification-by-kmer pipeline")
    args = parser.parse_args()

    jobtag = camisim.gen_prefix(args.readsGB)
    job_ids = {}
    log = get_log(jobtag, lvl=logging.DEBUG)

    ###### Fetch required ids. ######
    id_fp = os.path.join(working_dir, ncbi_id_file)
    if not os.path.isfile(id_fp):
        raise_w_log(log, IOError, "Expected id_list.txt not found -> "+id_fp)

    ncbi_ids = [entry.strip() for line in open(id_fp, "r").readlines() for entry in line.split()]
    log.info(f"Found {len(ncbi_ids)} ids in file")

    expected_output = [f"data/simulated_data/input_genomes/{x.replace('.', '_')}.fa" for x in ncbi_ids+["combined"]]
    if all([os.path.isfile(x) for x in expected_output]):
        log.info(f"Input genomes found - skipping load")
    else: 
        call_fetch_fastas = ncbifetch.generate_syscall(
            working_dir=working_dir, 
            ids=ncbi_ids,
            outdir=os.path.join(working_dir, 'data/simulated_data/input_genomes')
            )
        qsub_args = ncbifetch.gen_qsub_args(working_dir=working_dir)
        runid_fetch = submit2(command=call_fetch_fastas, test=False, **qsub_args)
        job_ids["runid_fetch"] = runid_fetch
        log.info("fetch jobid: " + runid_fetch)
    
    ###### Antismash input genomes. ######
    antismash_out = f"data/simulated_data/antismash/input_genomes/"
    antismash_html = os.path.join(antismash_out, "index.html")
    if os.path.isfile(antismash_html):
        log.info(f"Antismash for (input_genomes) exists - skipping...")
    else:
        combined_input_genomes = "data/simulated_data/input_genomes/combined.fa"
        antismash.preflight(fastafile=combined_input_genomes, outdir=antismash_out)
        call_antismash = antismash.generate_syscall(combined_input_genomes, antismash_out)
        antismash_qsub = antismash.gen_qsub_args(working_dir, dependency=job_ids.get("runid_fetch", []))
        runid_antismash_input = submit2(
            command = call_antismash,
            test=False,
            **antismash_qsub)
        job_ids["runid_antismash_input"] = runid_antismash_input
        log.info(f"Running Antismash Input Genomes jobid: {runid_antismash_input}")


    ############ Run CAMISIM ############
    if camisim.output_exists(args.readsGB):
        log.info(f"Camisim output for {args.readsGB}GB found - skipping...")
    else:
        log.info("Running Camisim")
        camisim.preflight(reads_GB=args.readsGB)
        camisim_qsub_kwargs = camisim.gen_qsub_args(job_tag = jobtag, dependency=job_ids.get("runid_fetch", []))
        runid_camisim = submit2(
            command=camisim.generate_syscall(reads_GB=args.readsGB),
            **camisim_qsub_kwargs
        )
        job_ids["runid_camisim"] = runid_camisim
        log.info(f"Running Camisim id {runid_camisim}")
    
    ############ Run QC pipeline  ############
    # TODO: Ask Trine


    ############ Run mapping quantification pipeline  ############
    # Dependent on antismash input and qc pipeline (for now camisim as standin)

    if not args.runQuantifierMap:
        log.info("Not running runQuantifierMap pipeline, set --runQuantifierFlag flag")
    else:
        readsfile = f"data/simulated_data/camisim/{jobtag}/sample_0/reads/anonymous_reads.fq.gz"
        reference = "data/simulated_data/antismash/input_genomes/combined_bgc.fa"
        outdir    = f"data/simulated_data/quantification_map/{jobtag}/"
        #Preppring and running job
        if QuantifierMap.is_success(outdir):
            log.info("Quantification-by-mapping Already run - skipping...")
        else:
            dependencies = [job_ids[x] for x in ("runid_camisim","runid_antismash_input") if x in job_ids]
            api = QuantifierMap(reads=readsfile, reference=reference, output_dir=outdir, minMapQ=30, log=log)
            api.preflight(check_input=False)
            api.set_qsub_args(jobtag=jobtag, dependency=dependencies)
            api.generate_syscall() #not needed as we run the default
            api.add_to_que(test=False)
            job_ids["QuantifierMap"] = api.job_id

    ############ Run KMER quantification pipeline (Simple) ############
    if not args.runQuantifierKmer:
        log.info("Not running runQuantifierKmer pipeline, set --runQuantifierKmer flag")
    else:
        readsfile   = f"data/simulated_data/camisim/{jobtag}/sample_0/reads/anonymous_reads.fq.gz"
        catalogue   = f"data/simulated_data/quantification_kmer/bgc_catalogues"
        outdir      = f"data/simulated_data/quantification_kmer/{jobtag}/"

        if QuantifierKmer.is_success(outdir):
            log.info("Quantification-by-kmer already run - skipping ...")
        else:
            if not os.path.isdir(catalogue):
                raise_w_log(log, IOError, "Catalogue files not found - must exists at init - check scripts/kmer_gen_catalogue for simple catalogue")
            api = QuantifierKmer(reads=readsfile, fp_catalogue=catalogue, output_dir=outdir, kmer_size = 15, log=log) #This should be in a config somewhere.
            api.preflight(check_input=False) 
            api.set_qsub_args(jobtag=jobtag)
            api.generate_syscall() #not needed as we run the default
            api.add_to_que()
            job_ids["QuantifierKmer"] = api.job_id

    ############ Run Assembly pipeline  ############ #TODO: Should we perhaps wrap steps in functions?
    if not args.runAssembly:
        log.info("Not running assembly pipeline, set --runAssembly flag")
    else:
        dir_p = f"data/simulated_data/camisim/{jobtag}/sample_0/reads"
        # r1_fn = ['Genome11.fq.gz','Genome21.fq.gz', 'Genome31.fq.gz','Genome41.fq.gz','Genome51.fq.gz']
        # r2_fn = ['Genome11.fq.gz','Genome21.fq.gz', 'Genome31.fq.gz','Genome41.fq.gz','Genome51.fq.gz']
        # reads1 = sorted([os.path.join(dir_p,x) for x in r1_fn])
        # reads2 = sorted([os.path.join(dir_p,x)  for x in r2_fn])
        reads_interleaved_fp = os.path.join(dir_p, "anonymous_reads.fq.gz" )
        assembly_outdir = f"data/simulated_data/assembly/{jobtag}"    
        #Check if already run:
        run_assembly = True
        if os.path.isdir(assembly_outdir):
            if Assembler.is_success(assembly_outdir):
                log.info("Assembly already run and is success - skipping")
                run_assembly=False
            else:
                log.info("Resetting outdir and rerunning assembly - previous run appears to have failed.")
                shutil.rmtree(assembly_outdir)
                run_assembly=True

        if run_assembly:
            api_Assembler = Assembler(reads_interleaved=reads_interleaved_fp, output_dir=assembly_outdir, log=log)
            api_Assembler.preflight(check_input=False) #Cant check input on dependency jobs as input might not exist at this point.
            api_Assembler.set_qsub_args(jobtag=jobtag, dependency = job_ids.get("runid_camisim",[]))
            api_Assembler.generate_syscall() #not needed as we run the default
            api_Assembler.add_to_que(test=False)
            runid_assembler = api_Assembler.job_id
            job_ids["runid_assembler"] = runid_assembler
            
        ############ Run Antismash on assembly  ############
        #IO
        antismash_out = f"data/simulated_data/antismash/{jobtag}/"
        antismash_html = os.path.join(antismash_out, "index.html")

        assembly_fp = os.path.join(assembly_outdir, "scaffolds.fasta")

        if os.path.isfile(antismash_html):
            log.info(f"Antismash for ({args.readsGB}) exists - skipping...")
        else:
            antismash.preflight(fastafile=assembly_fp, outdir=antismash_out)
            call_antismash = antismash.generate_syscall(assembly_fp, antismash_out)
            antismash_qsub = antismash.gen_qsub_args(working_dir, dependency=job_ids.get("runid_assembler",[]))
            runid_antismash_assembly = submit2(
                command = call_antismash,
                test=False,
                **antismash_qsub)
            log.info(f"Running Antismash assembly id {runid_antismash_assembly}")
    ########################
    if not job_ids:
        log.warning("No jobs were started for ->"+jobtag)
    log.info("Runid map (also stdout)"+str(job_ids))
    sys.stdout.write(str(job_ids)+"\n")

            
