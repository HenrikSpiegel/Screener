import os
import pathlib
import logging
from typing import Union
from pathlib import Path

from scripts.qsub_base import Base

class Antismash(Base):
    def __init__(self, fastafile: Union[str,Path], outdir: Union[str,Path],log: logging.Logger=None) -> None:
        if log:
            self.add_external_log(log)

        if not isinstance(fastafile, Path):
            fastafile = Path(fastafile)
        self.fatafile = fastafile
        if not isinstance(outdir, Path):
            outdir = Path(outdir)
        self.outdir = outdir

    qsub_requirements = dict(
        modules = "tools anaconda3/2021.11 antismash/6.1.1",
        runtime = 60,
        cores = 30,
        ram = 100,
        )

    def preflight(self, check_input=False) -> None:
        if check_input:
            if not self.fastafile.is_file():
                msg = f"Missing inputfile -> {self.fastafile}"
                self.log.error(msg)
                raise IOError(msg)
        self.outdir.mkdir(parents=True, exist_ok=True)

    def generate_syscall(self) -> None:
        # sets mem / cpu based on default qsub args.
        syscall = f"""\
antismash --output-dir {self.outdir} \
--taxon bacteria \
--cpus {self.qsub_requirements["cores"]-1} \
--genefinding-tool prodigal \
{self.fatafile}
#Pulls the output bgcs into 1 file.
python -m scripts.antismash_as_fasta -i {self.outdir}

touch {self.outdir / "success"}
"""
        self._syscall=syscall

    def successful(self):
        success_file = self.outdir / "success"
        return success_file.exists()

    @staticmethod
    def is_success(output_dir) -> str:
        success_file = Path(output_dir) / "success"
        return success_file.exists()