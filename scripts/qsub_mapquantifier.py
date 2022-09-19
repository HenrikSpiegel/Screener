import glob
import logging
import os
import pathlib
import re
import sys
from typing import List

from scripts.functions import submit2
from scripts.qsub_base import Base

class QuantifierMap(Base):
    def __init__(self, reads:str, reference:str, output_dir: str,log: logging.Logger=None) -> None:
        if log:
            self._log=log
        self.reads = self.as_abspath(reads)
        self.reference = self.as_abspath(reference)
        self.output_dir = self.as_abspath(output_dir)


    qsub_requirements = dict(
        modules = "tools anaconda3/2021.05 minimap2/2.6 samtools/1.14 cmseq/1.0.4",
        runtime = 360,
        cores = 38,
        ram=180,
        )

    def preflight(self, check_input=False) -> None:
        
        if check_input:
            if not os.path.isfile(self.reference):
                self.log.error("Reference file not found -> "+self.reference)
                raise IOError("Reference file not found -> "+self.reference)
            if not os.path.isfile(self.reads):
                self.log.error("Reads file not found -> "+self.reads)
                raise IOError("Reads file not found -> "+self.Reads)
            self.log.info("Found all input files")


        pathlib.Path(os.path.dirname(self.output_dir)).mkdir(parents=True, exist_ok=True)

    def generate_syscall(self) -> None:
        # sets mem / cpu based on default qsub args.

        # if not self.reference.endswith(".gbk"):
        syscall=f"""\
set -e        
minimap2 -t {self.qsub_args['cores']-2} -a {self.reference} {self.reads} > {os.path.join(self.output_dir, "aln.sam")} 
#Sort and index the resulting .sam file from the mapping.
samtools sort --threads {self.qsub_args['cores']-2} --write-index -o {os.path.join(self.output_dir, "sorted.aln.bam")} {os.path.join(self.output_dir, "aln.sam")}
samtools mpileup -a {os.path.join(self.output_dir, "sorted.aln.bam")} | awk '{{print $1"\t"$2"\t"$4}}' > {os.path.join(self.output_dir, "coverage.tsv")} 
#Get summation.
python3 -m cmseq.breadth_depth -f {os.path.join(self.output_dir, "sorted.aln.bam")} > {os.path.join(self.output_dir, "cmseq_summation.tsv")}
"""
        self._syscall = syscall
        return

#         self.log.info("reformatting reference .gbk to .fa before launching minimap2")
#         ref_fasta = os.path.basename(self.reference).replace(".gbk", ".fa")
#         ref_fasta = os.path.join(self.output_dir, ref_fasta)
#         syscall=f"""\
# python -m scripts.gbk2fa --gbk {self.reference} --fa {ref_fasta}
# minimap2 -t {self.qsub_args['cores']-2} -a {ref_fasta} {self.reads} > {os.path.join(self.output_dir, "aln.sam")} 
# #Sort and index the resulting .sam file from the mapping.
# samtools sort --threads {self.qsub_args['cores']-2} --write-index -o {os.path.join(self.output_dir, "sorted.aln.bam")} {os.path.join(self.output_dir, "aln.sam")}
# #Get summation.
# python3 -m cmseq.breadth_depth {os.path.join(self.output_dir, "sorted.aln.bam")} > summation.txt
# """
#         call_convert=f"python -m scripts.gbk2fa --gbk {self.reference} --fa {ref_fasta}"
#         syscall = "\n".join([call_convert, syscall])
#         self._syscall = syscall
        


    @staticmethod
    def is_success(output_dir) -> str:
        return os.path.isfile(os.path.join(output_dir, "cmseq_summation.tsv"))


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
    readsfile = "data/simulated_data/camisim/0_1GB/sample_0/reads/anonymous_reads.fq.gz"
    reference = "data/simulated_data/antismash/input_genomes/combined_bgc.fa"
    outdir    = "data/simulated_data/quantification/0_1GB/"
    #Preppring and running job
    api = QuantifierMap(reads=readsfile, reference=reference, output_dir=outdir)
    api.preflight(check_input=True)
    api.set_qsub_args(jobtag="test")
    api.generate_syscall() #not needed as we run the default
    api.add_to_que(test=False)

    print(api.job_id)