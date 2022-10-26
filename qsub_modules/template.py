#!/usr/bin/env python3
import logging
import os

default_qsub_args = dict(
    directory = "/home/projects/dtu_00009/people/henspi/git/AntibioticaScreening/project",
    modules = "tools, antismash/6.1.1",
    runtime = 120,
    cores = 30,
    ram=100,
    group="dtu_00009",
    jobname=os.path.basename(__file__),
    output = os.path.join("/home/projects/dtu_00009/people/henspi/git/AntibioticaScreening/project", "logs", os.path.basename(__file__), "_stdout"),
    error = os.path.join("/home/projects/dtu_00009/people/henspi/git/AntibioticaScreening/project", "logs", os.path.basename(__file__), "_stderr")
)

def prepare_run(log: logging.Logger = None, **kwargs):
    """Setup relevant files needed for the script to run.
    """
    ...

def main(log: logging.Logger = None, **kwargs) -> None:
    """perform the task.
    If the task is not 
    """
    ...

def _is_success(**kwargs) -> bool:
    """Check for whether the run has succesfully run.
    """
    ...

if __name__ == "__main__":
    prepare_run()
    main()