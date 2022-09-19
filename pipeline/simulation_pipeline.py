
import os, sys, time
import pathlib

from pipeline.general_functions import submit2
import scripts.qsub_camisim as camisim
import scripts.qsub_ncbiFetch as ncbifetch
import scripts.qsub_antismash as antismash

working_dir = "/home/projects/dtu_00009/people/henspi/git/AntibioticaScreening/project"
ncbi_id_file = "config/ids_simulation_genomes.txt"

camisim_dataset_GB_range = [1, 0.1, 0.01, 0.001, 0.0001, 0.00001]
camisim_dataset_GB_range = [0.1]


if __name__ == "__main__":
    
    ###### Fetch required ids. ######
    id_fp = os.path.join(working_dir, ncbi_id_file)
    if not os.path.isfile(id_fp):
        raise IOError("Expected id_list.txt not found -> "+id_fp)

    ncbi_ids = [entry.strip() for line in open(id_fp, "r").readlines() for entry in line.split()]
    print(f"Found {len(ncbi_ids)} ids in file", file=sys.stderr)
    call_fetch_fastas = ncbifetch.generate_syscall(
        working_dir=working_dir, 
        ids=ncbi_ids,
        outdir=os.path.join(working_dir, 'data/simulated_data/input_genomes')
        )
    qsub_args = ncbifetch.gen_qsub_args(working_dir=working_dir)
    runid_fetch = submit2(command=call_fetch_fastas, test=False, **qsub_args)
    print("fetch jobid: " + runid_fetch, file=sys.stderr)
    time.sleep(1)

    # TODO: Consider doing a direct check if finished instead of using que:
    # check if fastas script ran correctly
    # We could check if combined exists and has the correct size.
#     combined_fasta = 

    ###### Run antismash on the full genomes. ######
    fastafile = os.path.join(working_dir, "data/simulated_data/input_genomes/combined.fasta")
    antismash_out = os.path.join(working_dir, "data/simulated_data/input_genomes/antismash")
    antismash.preflight(fastafile=fastafile, outdir=antismash_out)
    antismash_runid_input = submit2(
        command = antismash.generate_syscall(fastafile, antismash_out),
        test=False,
        **antismash.gen_qsub_args(working_dir=working_dir, dependency=[runid_fetch])
    )
    print("antismash id (input genomes): "+antismash_runid_input, file=sys.stderr)
    
    ###### Generate synthetic metagenomic sample ######
    print("Starting CAMISIM jobs", file=sys.stderr)
    camisim_jobs = {}
    for reads_gb in camisim_dataset_GB_range:
        if camisim.output_exists(reads_gb):
            print(f"Output for {reads_gb}GB found - skipping...", file=sys.stderr)
            continue
        print(f"Starting job for CAMISIM: {reads_gb}GB", file=sys.stderr)
        camisim.preflight(reads_GB=reads_gb)
        jobtag = camisim.gen_prefix(reads_gb)
        camisim_qsub_kwargs = camisim.gen_qsub_args(job_tag = jobtag, dependency=[runid_fetch])
        runid = submit2(
            command=camisim.generate_syscall(reads_GB=reads_gb),
            **camisim_qsub_kwargs
        )
        camisim_jobs[camisim.gen_prefix(reads_gb)]=runid
    print(camisim_jobs, file=sys.stderr)

    ###### Run assembly on camisim jobs. ######

    # TODO: Setup assembly pipeline!

    ###### Run antismash on camisim re-assembled contigs. ######
    # TODO: For now we are simply using the GSA contigs as standin.
    # Check which assemblies needs antismashing
    antismash_jobs = {}
    for gb in camisim_dataset_GB_range:
        named_gb = camisim.gen_prefix(reads_gb)
        antismash_out = f"data/simulated_data/antismash/{named_gb}/"
        antismash_html = os.path.join(antismash_out, "index.html")
        if os.path.isfile(antismash_html):
            print(f"Antismash for ({gb}) exists - skipping...", file=sys.stderr)
            continue

        assembly_fp = f"/data/simulated_data/camisim/{named_gb}/*_sample_0/contigs/gsa.fasta.gz"
        call_unzip = f"gunzip {assembly_fp}"
        call_zip = f"gzip {assembly_fp}"

        antismash.preflight(fastafile=assembly_fp, outdir=antismash_out)

        call_antismash = antismash.generate_syscall(assembly_fp, antismash_out)
        combined_call = "\n".join([call_unzip, call_antismash, call_zip])
        antismash_qsub = antismash.gen_qsub_args(working_dir)
        
        camisim_dep = camisim_jobs.get(named_gb, []) #TODO should be assembler dependency
        if camisim_dep:
            print(f"{named_gb} in que - setting dependency", file=sys.stderr)
            dependency=camisim_jobs[named_gb]
        antismash_runid_input = submit2(
            command = combined_call,
            test=False,
            **antismash_qsub)





