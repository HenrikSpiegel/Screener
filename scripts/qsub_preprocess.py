from scripts.qsub_base import Base
import os
import pathlib
import glob
import logging
from typing import List, Union
from pathlib import Path


class Preprocessor(Base):
    def __init__(self, reads_interleaved:List[Path], outdir:Path, log: logging.Logger=None) -> None:
        """
        Note it is assummed that the grandparent +  file combination are unique.
        sample_0/reads/reads.fa.gz -> outdir/sample_0/trimmed.reads.fa.gz
        sample_1/reads/reads.fa.gz -> outdir/sample_1/trimmed.reads.fa.gz
        """
        if log:
            self.add_external_log(log)

        if not isinstance(reads_interleaved, list):
            reads_interleaved = [reads_interleaved]
        self.reads_files = [Path(x) for x in reads_interleaved]
        self.outdir = Path(outdir)

        #Here we assume that the grandparent is sample_id ie sample_0/reads/readsfile however we should perhaps redo
        self.fp_outs_interleaved = [self.outdir / infile.parent.parent.name / Path("trimmed."+infile.name) for infile in self.reads_files]
        self.fp_outs_singles     = [self.outdir / infile.parent.parent.name / Path("trimmed.single"+infile.name) for infile in self.reads_files]
                
    qsub_requirements = dict(
        modules = "tools sickle/20150314",
        runtime = 360,
        cores = 20,
        ram=80,
        )
    
    def preflight(self, check_input=False) -> None:
        if check_input:
            if not os.path.isfile(self.reads_file):
                err_msg = f"Reads file not found -> {self.reads_file}"
                self.log.error(err_msg)
                raise FileNotFoundError(err_msg)
        for outs in self.fp_outs_interleaved:
            outs.parent.mkdir(parents=True, exist_ok=True)


    def generate_syscall(self) -> None:
        # sets mem / cpu based on default qsub args.

        calls = [f"sickle pe --gzip-output -l 100 --pe-combo {infile} -t sanger -m {out_inter} -s {out_single}" for infile, out_inter, out_single in zip(self.reads_files, self.fp_outs_interleaved, self.fp_outs_singles)]
        syscall = "\n".join(calls)
        syscall += f"""
touch {self.outdir / "success"}
        """
        self._syscall = syscall
        return
    
    def successful(self):
        success_file = self.outdir / "success"
        return success_file.exists()


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