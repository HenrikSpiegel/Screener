import os
import pathlib
import logging
from typing import List

from scripts.qsub_base import Base

class NCBIFetch(Base):
    def __init__(self, id_list:List[str], outdir: str,log: logging.Logger=None) -> None:
        if log:
            self.add_external_log(log)

        self.ids = id_list
        self.outdir = outdir

    qsub_requirements = dict(
        modules = "tools anaconda3/2020.07",
        runtime = 20,
        cores = 1,
        ram = 10
        )

    def preflight(self) -> None:
        pathlib.Path(self.outdir).mkdir(parents=True, exist_ok=True)

    def generate_syscall(self) -> None:
        # sets mem / cpu based on default qsub args.
        syscall = f"""\
python {"scripts/ncbi_fetch.py"} \
--ncbi_ids {" ".join(self.ids)} \
--outdir {self.outdir}\
"""
        self._syscall=syscall

    @staticmethod
    def is_success(output_dir) -> str:
        raise NotImplementedError