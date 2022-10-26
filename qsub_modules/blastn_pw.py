
import time
from qsub_modules.base import Base, QsubStatus
from Bio import SeqIO
from typing import Union, List
from pathlib import Path

class PairwiseBlast(Base):
    def __init__(self, fasta_file: Path, output_dir: Path):
        self.fasta_file = Path(fasta_file)

        self.output_dir = Path(output_dir)
        self.output_combined = Path(output_dir)/"combined_blast_results.tsv"

        self.success_file = self.output_dir/'success'

        self.output_columns = "qaccver saccver pident length qcovhsp mismatch evalue bitscore"

    qsub_requirements = dict(
        modules = "tools anaconda3/2021.05 perl/5.30.2 ncbi-blast/2.12.0+",
        runtime = 30,
        cores = 10,
        ram=20,
        )

    def preflight(self, check_input=True):
        self.output_dir.mkdir(parents=True, exist_ok=True)

        #prepare outputfile
        self.output_combined.unlink(missing_ok=True)
        self.output_combined.write_text(self.output_columns.replace(" ", "\t"))

        self.success_file.unlink(missing_ok=True)

        if check_input:
            assert(self.fasta_file.exists())

    # def prepare_pairwise_files(self):
    #     self.entries=[]
    #     for file in self.fasta_file:
    #         for record in SeqIO.parse(file, "fasta"):
    #             outfile = self.fasta_dir/(record.name+".fa")
    #             self.entries.append(outfile)
    #             if outfile.is_file():
    #                 continue
    #             self.fasta_file.append(outfile)
    #             entry = f">{record.description}\n{record.seq}"
    #             outfile.write_text(entry)
    #     self.log.info(f"Individual fastafiles added to {self.fasta_dir}")

    def generate_syscall(self):

        call = f"""\
blastn -query {self.fasta_file} -subject {self.fasta_file} -task dc-megablast -outfmt "6 {self.output_columns}" -max_hsps 1 -evalue 5000 >> {self.output_combined}

touch {self.success_file}
"""
        self._syscall=call

if __name__ == "__main__":
    
    pb = PairwiseBlast(
        "data/simulated_data/antismash/input_genomes/combined_bgc.fa", 
        "data/simulated_data/pairwise")
    pb.preflight()
    pb.add_to_que()

    pb.wait_for_finish()

    
