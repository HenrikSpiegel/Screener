import glob
import logging
import os
import pathlib
import re
import sys
from typing import List

from scripts.functions import submit2
from scripts.qsub_base import Base

class Assembler(Base):
    def __init__(self, reads_interleaved: str, output_dir: str="data/simulated_data/assembly/default", log: logging.Logger=None) -> None:
        if log:
            self._log=log

        self.output_dir = output_dir
        self.reads_interleaved = reads_interleaved

    qsub_requirements = dict(
        modules = "tools anaconda3/2021.05 spades/3.15.2",    #"tools megahit/1.2.9",
        runtime = 360,
        cores = 38,
        ram=180,
        )

    def preflight(self, check_input=False) -> None:
        self.log.debug("running preflight")
        if check_input: #cant be used for qued jobs.
            self.log.info("checking input")
            assert(os.path.isfile(self.reads_interleaved))
            # all_input = self.reads1+self.reads2
            # files_missing = [x for x in all_input if not os.path.isfile(x) ]
            # if files_missing:
            #     self.log.error(f"Input files not found: {files_missing}")
            #     raise ValueError(f"Input files not found: {files_missing}")
        #could check if output dir is filled.
        if os.path.isdir(self.output_dir):
            if self.is_success(self.output_dir):
                self.log.info("Previous run appears succesful.")
            raise ValueError(f"Output directory already exist - cannot overwrite -> {self.output_dir}")
        pathlib.Path(os.path.dirname(self.output_dir)).mkdir(parents=True, exist_ok=True)

    def generate_syscall(self) -> str:
        # sets mem / cpu based on default qsub args.
#         syscall = f"""\
# megahit \
# -1 {",".join(self.reads1)} \
# -2 {",".join(self.reads2)} \
# --memory {self.qsub_args["ram"]-5} \
# --num-cpu-threads {self.qsub_args["cores"]-1} \
# --out-dir {self.output_dir}\
# """
        syscall = f"""\
python /services/tools/spades/3.15.2/bin/metaspades.py \
--12 {self.reads_interleaved} \
-o {self.output_dir} \
--threads {self.qsub_args["cores"]-1} \
--memory {self.qsub_args["ram"]-5} \
-k 27,47,67,87,107,127 \
"""
        self._syscall = syscall
    
    @staticmethod
    def is_success(output_dir) -> str:
        return os.path.isfile(os.path.join(output_dir, "scaffolds.fasta"))


if __name__ == "__main__":
    #IO
    gb = 0.1

    # fussy_dir = f"data/simulated_data/camisim/{Base.gen_prefix(gb)}/sample_0/reads"
    # dir_p = glob.glob(fussy_dir)[0]
    # files = os.listdir(dir_p)
    # reads1 = sorted([os.path.join(dir_p,x) for x in files if re.search("Genome\d1.fq.gz", x)])
    # reads2 = sorted([os.path.join(dir_p,x)  for x in files if re.search("Genome\d2.fq.gz", x)])
    reads = f"data/simulated_data/camisim/{Base.gen_prefix(gb)}/sample_0/reads/anonymous_reads.fq.gz"
    
    outdir = f"data/simulated_data/assembly/{Base.gen_prefix(gb)}"
    
    #Preppring and running job
    api = Assembler(reads_interleaved=reads, output_dir=outdir)
    api.preflight(check_input=True)
    api.set_qsub_args(jobtag=Base.gen_prefix(gb))
    api.generate_syscall() #not needed as we run the default
    api.add_to_que(test=False)

    print(api.job_id)