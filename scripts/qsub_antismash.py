#!/usr/bin/env python3
import logging
import os, sys
import pathlib

from scripts.functions import submit2

def gen_qsub_args(working_dir, **kwargs):
    scriptname = os.path.basename(__file__).split(".")[0]
    qsub_args = dict(
        directory = working_dir,
        modules = "tools antismash/6.1.1",
        runtime = 60,
        cores = 30,
        ram = 100,
        group = "dtu_00009",
        jobname=scriptname,
        output = os.path.join(working_dir, "logs", scriptname+ "_stdout"),
        error = os.path.join(working_dir, "logs", scriptname+ "_stderr")
    )
    qsub_args.update(kwargs)
    return qsub_args

def generate_syscall(fastafile: str, outdir:str=None, log: logging.Logger = None) -> None:
    syscall = f"""\
antismash --output-dir {outdir} \
--taxon bacteria \
--cpus 30 \
--genefinding-tool prodigal \
{fastafile}
#Pulls the output bgcs into 1 file.
python -m scripts.antismash_as_fasta -i {outdir}
"""
    if log: log.info("running:\n"+syscall)
    return syscall

def preflight(fastafile: str, outdir: str,  log: logging.Logger = None):
    """Setup relevant files needed for the script to run.
    """
    pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)
    #if os.path.isfile(fastafile):
    #    if log: log.debug("Input file detected")
    #    return
    #if log: log.error("Input file not found -> "+fastafile)
    #raise IOError("Input file not found -> "+fastafile)




def _is_success(**kwargs) -> bool:
    """Check for whether the run has succesfully run.
    """
    ...

if __name__ == "__main__":
    import argparse
    ## Front matter - handle input parsing
    parser = argparse.ArgumentParser()
    parser.add_argument('--fastafile', required=True, help="path to fastafile containing 1 or more entries.")
    parser.add_argument('--outdir', default=None, help="target for outdir, default subdir of fastafile")
    args = parser.parse_args()
    working_dir = "/home/projects/dtu_00009/people/henspi/git/AntibioticaScreening/project"

    if args.outdir is None:
        args.outdir = args.fastafile.split(".")[0] + "_antismash"
        print("outdir set to -> "+args.outdir, file=sys.stderr)

    preflight(fastafile=args.fastafile, outdir=args.outdir)
    runid = submit2(
        command = generate_syscall(args.fastafile, args.outdir),
        test=False,
        **gen_qsub_args(working_dir=working_dir)
    )
    print("antismash id: "+runid, file=sys.stderr)