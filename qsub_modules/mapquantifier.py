import glob
import logging
import os
import pathlib
from pathlib import Path
import re
import sys
from typing import List, Union

from qsub_modules.functions import submit2
from qsub_modules.base import Base

class QuantifierMap(Base):
    def __init__(self, reads:Union[str, List[str]], reference:str, output_dir: str, minMapQ: int=30, minBaseQ:int=20, log: logging.Logger=None) -> None:
        if log:
            self.add_external_log(log)
        
        if isinstance(reads, str):
            reads = [reads]
        self.reads = [Path(x) for x in reads]
        self.reference  = Path(reference)
        self.output_dir = Path(output_dir)
        self.minMapQ    = minMapQ
        self.minBaseQ   = minBaseQ


    qsub_requirements = dict(
        modules = "tools anaconda3/2021.05 minimap2/2.6 samtools/1.14 cmseq/1.0.4",
        runtime = 360,
        cores = 38,
        ram=180,
        )

    def preflight(self, check_input=False) -> None:
        
        if check_input:
            if not self.reference.is_file():
                self.log.error("Reference file not found -> "+self.reference)
                raise IOError("Reference file not found -> "+self.reference)
            if not all(x.is_file() for x in self.reads):
                self.log.error("Reads file not found -> "+self.reads)
                raise IOError("Reads file not found -> "+self.Reads)
            self.log.info("Found all input files")

        self.output_dir.mkdir(parents=True, exist_ok=True)

    def generate_syscall(self) -> None:
        syscall=f"""\
minimap2 -t {self.qsub_args['cores']-2} -a {self.reference} {" ".join(x.__str__() for x in self.reads)} > {self.output_dir / "aln.sam"} 

#Sort and index the resulting .sam file from the mapping.
samtools sort --threads {self.qsub_args['cores']-2} --write-index -o {self.output_dir / "sorted.aln.bam"} {self.output_dir / "aln.sam"}

#Get by position coverage
samtools mpileup --min-MQ {self.minMapQ} --min-BQ {self.minBaseQ} -a {self.output_dir / "sorted.aln.bam"} | awk '{{print $1"\t"$2"\t"$4}}' > {self.output_dir / "coverage.tsv"}
samtools mpileup --min-MQ 0 --min-BQ 0 -a {self.output_dir / "sorted.aln.bam"} | awk '{{print $1"\t"$2"\t"$4}}' > {self.output_dir / "coverage_raw.tsv"}  
#Get summation.
python3 -m cmseq.breadth_depth --minqual {self.minBaseQ} -f {self.output_dir / "sorted.aln.bam"} > {self.output_dir / "cmseq_summation.tsv"}
python3 -m cmseq.breadth_depth --minqual 0 -f {self.output_dir / "sorted.aln.bam"} > {self.output_dir / "cmseq_summation_raw.tsv"}
"""
        self._syscall = syscall
        return

    def successful(self):
        success_file = self.output_dir / "cmseq_summation.tsv"
        if not success_file.is_file():
            return False
        return success_file.read_bytes().count(b"\n") > 2

    @staticmethod
    def is_success(output_dir) -> str:
        return (Path(output_dir) / "cmseq_summation.tsv").is_file()


if __name__ == "__main__":
    # import argparse
    # ## Front matter - handle input parsing
    # parser = argparse.ArgumentParser()
    # parser.add_argument("--reads", required=True, help="file containing reads to be mapped")
    # parser.add_argument("--ref", required=True, help="reference db")
    # parser.add_argument("-o", required=True, help="output dir")
    # args = parser.parse_args()

    # readsfile = args.reads
    # referece  = args.ref
    # outdir    = args.o

    #IO
    readsfile = "data/simulated_data/camisim/0_01GB/sample_0/reads/anonymous_reads.fq.gz"
    reference = "data/simulated_data/antismash/input_genomes/combined_bgc.fa"
    outdir    = "data/simulated_data/quantification_map/0_1GB/"
    #Preppring and running job
    api = QuantifierMap(reads=readsfile, reference=reference, output_dir=outdir)
    api.preflight(check_input=True)
    api.set_qsub_args(jobtag="test")
    api.generate_syscall() #not needed as we run the default
    api.add_to_que(test=False)

    print(api.job_id)