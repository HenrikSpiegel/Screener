
from scripts.qsub_base import Base
from Bio import SeqIO
from typing import Union, List
from pathlib import Path

class PairwiseBlast(Base):
    def __init__(self, fasta_file: Union[Path, List[Path]], output_dir: Path):
        if isinstance(fasta_file, List):
            self.fasta_file = [Path(x) for x in fasta_file]
        else:
            self.fasta_file = [Path(fasta_file)]

        assert(all(file.exists() for file in self.fasta_file))

        self.output_dir = Path(output_dir)
        self.fasta_dir  = Path(output_dir)/"fastas"

        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.fasta_dir.mkdir(parents=True, exist_ok=True)

        self.fasta_files = []
        self.pairwise_files = []

    
    qsub_requirements = dict(
        modules = "tools anaconda3/2021.05 perl/5.30.2 ncbi-blast/2.12.0+",
        runtime = 30,
        cores = 10,
        ram=20,
        )

    def preflight(self):
        self.prepare_pairwise_files()
        self.combinations = [(self.entries[i], self.entries[i+1::]) for i in range(len(self.entries)-1)]
        calls = []
        for comb in self.combinations:
            f1 = comb[0]
            for f2 in comb[1]:
                out_fn = self.output_dir /"blast"/ (f1.stem+"-"+f2.stem)
                self.pairwise_files.append(out_fn)

            call = f"blastn -query {f1} -subject {f2} -outfmt 6 -max_hsps 1 > {out_fn}"
            calls.append(call)
        self.calls = calls

    def prepare_pairwise_files(self):
        self.entries=[]
        for file in self.fasta_file:
            for record in SeqIO.parse(file, "fasta"):
                outfile = self.fasta_dir/(record.name+".fa")
                self.entries.append(outfile)
                if outfile.is_file():
                    continue
                self.fasta_file.append(outfile)
                entry = f">{record.description}\n{record.seq}"
                outfile.write_text(entry)

    def generate_syscall(self):

        call = "\n".join(self.calls)
        call += f"\ntouch {self.output_dir/'success'}"
        self._syscall=call
