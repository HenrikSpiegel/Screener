
import shutil
import time
from qsub_modules.base import Base, QsubStatus
from Bio import SeqIO
from typing import Union, List
from pathlib import Path

class PairwiseBlast(Base):
    def __init__(self, fasta_file: Path, output_dir: Path, log=None, loglvl = "DEBUG"):

        if log:
            self.add_external_log(log)
        self.loglvl = loglvl

        self.fasta_file = Path(fasta_file)

        self.output_dir = Path(output_dir)
        self.output_combined = Path(output_dir)/"combined_blast_results.tsv"

        self.success_file = self.output_dir/'success'

        self.database_path  = self.output_dir / "blastDB" / self.fasta_file.stem
        self.output_columns = "qaccver saccver qstart qend sstart send pident length qcovs qcovhsp mismatch evalue bitscore"

    qsub_requirements = dict(
        modules = "tools anaconda3/2021.05 perl/5.30.2 ncbi-blast/2.12.0+",
        runtime = 360,
        cores = 38,
        ram=180,
        )

    def preflight(self, check_input=True):
        self.output_dir.mkdir(parents=True, exist_ok=True)

        #prepare outputfile
        self.output_combined.unlink(missing_ok=True)
        self.output_combined.write_text(self.output_columns.replace(" ", "\t")+"\n")

        #prepare for db:
        if self.database_path.parent.is_dir():
            self.log.info("Removing old DB")
            shutil.rmtree(self.database_path.parent)
        self.database_path.parent.mkdir(parents=True)

        if check_input:
            assert(self.fasta_file.exists())

    def generate_syscall(self):
        
        self._syscall = f"""\
makeblastdb -in {self.fasta_file} -parse_seqids -max_file_sz '4GB' -blastdb_version 5 -title "Database" -dbtype nucl -out {self.database_path}

blastn -query {self.fasta_file} -db {self.database_path} -task dc-megablast -outfmt "6 {self.output_columns}" -num_threads {self.qsub_args["cores"]-1} -subject_besthit >> {self.output_combined}

python scripts/symmetrise_blastn.py --blast {self.output_combined} --fasta {self.fasta_file}
"""

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fasta", required=True, type=Path, help="Multifasta file containing all partners for the pairwise blast")
    parser.add_argument("-o", required=True, type=Path, help="output dir")
    args = parser.parse_args()

    pb = PairwiseBlast(
        fasta_file=args.fasta,
        output_dir=args.o)
    pb.preflight()
    pb.add_to_que()


    
