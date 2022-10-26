import glob
import logging
import os
import pathlib
import re
import sys
from typing import List

from qsub_modules.functions import submit2
from qsub_modules.qsub_base import Base

class FastQC(Base):
    def __init__(self, reads:str, output_dir: str,log: logging.Logger=None) -> None:
        if log:
            self._log=log
        self.reads = reads
        self.output_dir = self.as_abspath(output_dir)

    qsub_requirements = dict(
        modules = "tools perl/5.30.2 jdk/18.0.1 fastqc/0.11.9",
        runtime = 30,
        cores = 10,
        ram=40,
        )

    def preflight(self, check_input=False) -> None:
        if check_input:
            if not os.path.isfile(self.reads):
                self.log.error("readfile file not found -> "+self.reads)
                raise IOError("readfile file not found -> "+self.reads)
        pathlib.Path(self.output_dir).mkdir(parents=True, exist_ok=True)

    def generate_syscall(self) -> None:
        # sets mem / cpu based on default qsub args.

        if self.reads.endswith(".gz"):
            infile = self.reads[:-3]
            syscall=f"""\
gunzip {self.reads}
fastqc -o {self.output_dir} -t {self.qsub_args['cores']-2} {infile}
gzip {infile}\
"""
        else:
            syscall=f"fastqc -o {self.output_dir} -t {self.qsub_args['cores']-2} {infile}"
        self._syscall=syscall

    @staticmethod
    def is_success(output_dir) -> str:
        raise NotImplementedError


if __name__ == "__main__":
    import argparse
    ## Front matter - handle input parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("--reads", required=True, help="file containing reads to be QC'ed")
    parser.add_argument("-o", required=True, help="output dir")
    args = parser.parse_args()

    readsfile = args.reads
    outdir    = args.o

    #Preppring and running job
    api = FastQC(reads=readsfile, output_dir=outdir)
    api.preflight(check_input=True)
    api.set_qsub_args(jobtag="M")
    api.generate_syscall() #not needed as we run the default
    api.add_to_que(test=False)

    print(api.job_id)