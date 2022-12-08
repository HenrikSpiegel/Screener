
import shutil
import time
from qsub_modules.base import Base, QsubStatus
from Bio import SeqIO
from typing import Union, List
from pathlib import Path

class MAGinator(Base):
    def __init__(self, MAGinator_dir:Path, MAGinator_wd: Path, refined_set_size=100,log=None, loglvl = "DEBUG"):

        if log:
            self.add_external_log(log)
        self.loglvl = loglvl
        self.MAGinator_dir=Path(MAGinator_dir) #="/home/projects/dtu_00009/people/henspi/git/MAGinator/"

        self.snake_pre = self.MAGinator_dir    / "maginator/workflow/prescreening_genes.Snakefile"
        self.snake_refine = self.MAGinator_dir / "maginator/workflow/signature_genes.Snakefile"
        self.ref_set_size = refined_set_size
        self.MAGinator_datadir =MAGinator_wd  #"/home/projects/dtu_00009/people/henspi/git/Screener/data/simulated_data/MAGinator"

        self.success_file = self.MAGinator_dir/ ".sucess_main"


    qsub_requirements = dict(
        modules = "tools anaconda3/2021.05 mamba-org/mamba/0.24.0",
        runtime = 12*60,
        cores = 38,
        ram=180,
        )

    def preflight(self, check_input=True):
        for file in (self.snake_pre, self.snake_refine, self.MAGinator_datadir):
            if not file.exists():
                raise FileNotFoundError(file.as_posix())
        #self.output_dir.mkdir(parents=True, exist_ok=True)


    def generate_syscall(self):
        
        self._syscall = f"""\
# part of config reads=reads.csv contigs=contigs.fasta vamb=clusters.tsv params=test/parameters.tab

snakemake --use-conda -s "{self.snake_pre}" --resources mem_gb={self.qsub_requirements["ram"]} --config wd="{self.MAGinator_datadir}" n_refined_genes={self.ref_set_size}\
 contigs="{self.MAGinator_datadir}/unused/contigs.fasta" vamb="{self.MAGinator_datadir}/unused/clusters.tsv" params="{self.MAGinator_datadir}/unused/parameters.tab"\
  --cores {self.qsub_requirements["cores"]} --printshellcmds format_conversion prescreening_genes

snakemake --use-conda -s "{self.snake_refine}" --resources mem_gb={self.qsub_requirements["ram"]} \
--config wd="{self.MAGinator_datadir}" n_refined_genes={self.ref_set_size} contigs="{self.MAGinator_datadir}/unused/contigs.fasta" vamb="{self.MAGinator_datadir}/unused/clusters.tsv" params="{self.MAGinator_datadir}/unused/parameters.tab"\
 --cores {self.qsub_requirements["cores"]} --printshellcmds --forcerun --until refinement

"""
