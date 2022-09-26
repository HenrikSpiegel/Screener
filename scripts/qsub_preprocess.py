from scripts.qsub_base import Base
import os
import pathlib
import glob
import logging
from typing import List, Union


class Preprocessor(Base):
    def __init__(self, reads_interleaved:str, outdir, log: logging.Logger=None) -> None:
        if log:
            self.add_external_log(log)

        self.reads_file = reads_interleaved
        self.outdir = outdir

        fname = os.path.basename(reads_interleaved)
        self.fp_out_interleaved = os.path.join(outdir, "trimmed."+fname)
        self.fp_out_singles     = os.path.join(outdir, "trimmed.single."+fname)        
        
    qsub_requirements = dict(
        modules = "tools sickle/20150314",
        runtime = 60,
        cores = 20,
        ram=80,
        )
    
    @property
    def log_setup(self):
        setup = super().log_setup
        setup.update({"level":logging.INFO})
        return setup

    def preflight(self, check_input=False) -> None:
        if check_input:
            if not os.path.isfile(self.reads_file):
                err_msg = f"Reads file not found -> {self.reads_file}"
                self.log.error(err_msg)
                raise FileNotFoundError(err_msg)
        pathlib.Path(self.outdir).mkdir(parents=True, exist_ok=True)


    def generate_syscall(self) -> None:
        # sets mem / cpu based on default qsub args.

        # if not self.reference.endswith(".gbk"):
        syscall=f"""
set -e
sickle pe --gzip-output -l 100 -c {self.reads_file} -t sanger -m {self.fp_out_interleaved} -s {self.fp_out_singles}
    """
        self._syscall = syscall
        return

    @staticmethod
    def is_success(outdir) -> str:
        possible_outputs = glob.glob(os.path.join(outdir, "trimmed.*"))
        return len(possible_outputs) == 2

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", required=True, help="file containing interleaved fq files")
    parser.add_argument("-o", required=True, help="output dir")
    args = parser.parse_args()

    api = Preprocessor(reads_interleaved=args.f, outdir=args.o)
    api.preflight(check_input=True)
    api.set_qsub_args(jobtag="test")
    api.generate_syscall() #not needed as we run the default
    #print(api.syscall)
    api.add_to_que(test=False)