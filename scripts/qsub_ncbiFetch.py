#!/usr/bin/env python3
import os, sys
import pathlib
from subprocess import call
from typing import List, Union

from scripts.functions import submit2

def gen_qsub_args(working_dir, **kwargs):
    scriptname = os.path.basename(__file__).split(".")[0]
    qsub_args = dict(
        directory = working_dir,
        modules = "tools anaconda3/2020.07",
        runtime = 20,
        cores = 1,
        ram = 10,
        group = "dtu_00009",
        jobname=scriptname,
        output = os.path.join(working_dir, "logs", scriptname+ "_stdout"),
        error = os.path.join(working_dir, "logs", scriptname+ "_stderr")
    )
    qsub_args.update(kwargs)
    return qsub_args

def generate_syscall(working_dir: str, ids: List[str], outdir: str = None):
    if outdir is None:
        outdir = os.path.join(working_dir, 'data/simulated_data/input_genomes')
        print("Outdir set to: "+outdir, file=sys.stderr)
    call_fetchncbi = f"""\
python {os.path.join(working_dir, "scripts/ncbi_fetch.py")} \
--ncbi_ids {" ".join(ids)} \
--outdir {outdir}\
"""
    return call_fetchncbi

def preflight(outdir: str):
    if not os.path.isdir(outdir):
        print("Creating outdir -> "+outdir, file=sys.stderr)
        pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)


if __name__ == "__main__":
    working_dir = "/home/projects/dtu_00009/people/henspi/git/AntibioticaScreening/project"
    id_list = ["NZ_LT906445.1","NZ_CP053893.1","NC_014328.1","NZ_CP020566.1","NZ_LT906470.1"]
    preflight(working_dir)
    qsub_kwargs = gen_qsub_args(working_dir=working_dir)
    runid_fetch = submit2(command=generate_syscall(working_dir=working_dir, ids=id_list), **qsub_kwargs)

    print("ncbifetch id: " + runid_fetch, file=sys.stderr)