
import shutil
import time
from qsub_modules.base import Base, QsubStatus
from Bio import SeqIO
from typing import Union, List
from pathlib import Path

class PairwiseFastANI(Base):
    def __init__(self, fasta_file: Path, output_dir: Path):
        self.fasta_file = Path(fasta_file)

        self.output_dir = Path(output_dir)
        self.output_combined = Path(output_dir)/"c.tsv"

        self.success_file = self.output_dir/'success'

        if not make_database:
            self.log.warning("Running the pairwise as query=file, subject=file which is slow for larger comparisons.")
        self.make_database = make_database
        self.database_path  = self.output_dir / "blastDB" / self.fasta_file.stem
        self.output_columns = "qaccver saccver qstart qend sstart send pident length qcovs qcovhsp mismatch evalue bitscore"

    qsub_requirements = dict(
        modules = "tools snakemake/7.14.0 fastani/1.33", #Not that snakemake and not python is required as some c++ libraries are silently needed.
        runtime = 120,
        cores = 30,
        ram=100,
        )

    def preflight(self, check_input=True):
        self.output_dir.mkdir(parents=True, exist_ok=True)

        #prepare outputfile
        self.output_combined.unlink(missing_ok=True)
        self.output_combined.write_text(self.output_columns.replace(" ", "\t")+"\n")

        #prepare for db:
        if self.make_database:
            if self.database_path.parent.is_dir():
                self.log.info("Removing old DB")
                shutil.rmtree(self.database_path.parent)
            self.database_path.parent.mkdir(parents=True)

        if check_input:
            assert(self.fasta_file.exists())

    def generate_syscall(self):
        
        if self.make_database:
            call = f"""\
makeblastdb -in {self.fasta_file} -parse_seqids -blastdb_version 5 -title "Database" -dbtype nucl -out {self.database_path}

blastn -query {self.fasta_file} -db {self.database_path} -task dc-megablast -outfmt "6 {self.output_columns}" -num_threads {self.qsub_args["cores"]-1} >> {self.output_combined}

touch {self.success_file}
"""
        else:
            call = f"""\
blastn -query {self.fasta_file} -subject {self.fasta_file} -task dc-megablast -outfmt "6 {self.output_columns}" -num_threads {self.qsub_args["cores"]-1} >> {self.output_combined}
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

    
