#!/usr/bin/env python3
import logging
import os, sys
import pathlib
from typing import List, Union

from qsub_modules.functions import submit2

def gen_qsub_args(working_dir:str=None, job_tag:str="", **kwargs):
    if working_dir is None:
        working_dir=os.getcwd()

    scriptname = os.path.basename(__file__).split(".")[0]
    scriptname = "_".join([scriptname, job_tag])
    qsub_args = dict(
        directory = working_dir,
        modules = "tools emboss/6.6.0",
        runtime = 120,
        cores = 10,
        ram = 50,
        group = "dtu_00009",
        jobname=scriptname,
        output = os.path.join(working_dir, "logs", scriptname+ "_stdout"),
        error = os.path.join(working_dir, "logs", scriptname+ "_stderr")
    )
    qsub_args.update(kwargs)
    return qsub_args

def generate_syscall(seq1_fp: str, seq2_fp: Union[List[str], str], outfile: str) -> None:
    if isinstance(seq2_fp, list): 
        seq2_fp = " ".join(seq2_fp)
    syscall = f"""\
needle \
-asequence {seq1_fp} \
-bsequence {seq2_fp} \
-outfile {outfile} \
-gapopen 10 \
-gapextend 0.5\
"""
    return syscall

def preflight(seq1_fp: str, seq2_fp: Union[List[str],str], outfile: str):
    """Setup relevant files needed for the script to run.
    """
    if not os.path.isfile(seq1_fp):
        raise IOError(f"seq1_fp not found -> " + seq1_fp)
    
    if not isinstance(seq2_fp, list):
        seq2_fp = [seq2_fp]
    for p in seq2_fp:
        if not os.path.isfile(p):
            raise IOError(f"seq2_fp not found -> " + seq1_fp)
    pathlib.Path(os.path.dirname(outfile)).mkdir(parents=True, exist_ok=True)


def _is_success(**kwargs) -> bool:
    """Check for whether the run has succesfully run.
    """
    ...

if __name__ == "__main__":
    import argparse
    ## Front matter - handle input parsing
    parser = argparse.ArgumentParser()
    parser.add_argument('--seq1', required=True, help="path to first file")
    parser.add_argument('--seq2', required=True, help="path to sencond file")
    parser.add_argument('--outfile', required=True, help="target for outputfile")
    args = parser.parse_args()
    working_dir = "/home/projects/dtu_00009/people/henspi/git/AntibioticaScreening/project"
   # print(generate_syscall(seq1_fp=args.seq1, seq2_fp=args.seq2, outfile=args.outfile))
    preflight(seq1_fp=args.seq1, seq2_fp=args.seq2, outfile=args.outfile)
    runid = submit2(
        command = generate_syscall(seq1_fp=args.seq1, seq2_fp=args.seq2, outfile=args.outfile),
        test=True,
        **gen_qsub_args(working_dir=working_dir)
    )
    print("needle id: "+runid, file=sys.stderr)