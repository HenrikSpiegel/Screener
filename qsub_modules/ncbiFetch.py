import os
import pathlib
import logging
from typing import List
from pathlib import Path

from qsub_modules.base import Base

class NCBIFetch(Base):
    def __init__(self, acc_ids:list[str], tax_ids:List[int], outdir: str,log: logging.Logger=None, loglvl = "DEBUG") -> None:
        if log:
            self.add_external_log(log)
        self.loglvl = loglvl
        self.acc_ids = acc_ids
        self.tax_ids = tax_ids
        self.outdir = Path(outdir)
        self.success_file = Path(outdir)/'success'

    qsub_requirements = dict(
        modules = "tools anaconda3/2020.07",
        runtime = 20,
        cores = 1,
        ram = 10
        )

    def preflight(self, check_input=True) -> None:
        pathlib.Path(self.outdir).mkdir(parents=True, exist_ok=True)

    def generate_syscall(self) -> None:
        # sets mem / cpu based on default qsub args.
        syscall = f"""\
python {"scripts/ncbi_fetch.py"} \
--acc_id {" ".join(self.acc_ids)} \
--tax_id {" ".join(self.tax_ids)} \
--outdir {self.outdir}\
"""
        self._syscall=syscall

class NCBIFetch_Assemblies(Base):
    def __init__(self, genera_names:List[str], genera_table_dir:Path, output_dir:Path, log: logging.Logger=None, loglvl = "DEBUG") -> None:
        if log:
            self.add_external_log(log)
        self.loglvl = loglvl

        self.genera_names = genera_names

        self.genera_table_dir   = Path(genera_table_dir)
        self.output_dir         = Path(output_dir)
        self.genomes_dir        = self.output_dir / 'genomes'
        self.taxdump_dir        = self.output_dir / "taxdump"
        self.taxdump            = self.taxdump_dir / 'ncbi_tax_dump.tar.gz'
        self.genera_table_fps   = [self.genera_table_dir / (genera.replace(" ","_") +".tsv") for genera in genera_names]

        self.assembly_dump_fp   = self.genera_table_dir / "assembly_summary.txt"

        self.success_file = self.output_dir / ".success_fetch"

    qsub_requirements = dict(
        modules = "tools anaconda3/2020.07",
        runtime = 60,
        cores = 2,
        ram = 20
        )

    def preflight(self, check_input=True) -> None:
        self.genera_table_dir.mkdir(parents=True, exist_ok=True)
        self.genomes_dir.mkdir(parents=True, exist_ok=True)
        self.taxdump_dir.mkdir(parents=True, exist_ok=True)

        self.genera_tables_exists = all(x.is_file() for x in self.genera_table_fps)
        


    def generate_syscall(self) -> None:
        # sets mem / cpu based on default qsub args.
        call_assembly_dump = f"""
wget --directory-prefix {self.genera_table_dir} https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt

#loading NCBI taxdump.
wget --directory-prefix {self.taxdump_dir} https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
mkdir {self.taxdump_dir}/NCBI
tar -xvf {self.taxdump_dir}/taxdump.tar.gz -C {self.taxdump_dir}/NCBI #For some fucking reason the archieve must be wrapped in a NCBI folder for CAMISIM to get going.
tar -czvf {self.taxdump} {self.taxdump_dir}/NCBI
rm -rf {self.taxdump_dir}/NCBI
rm {self.taxdump_dir}/taxdump.tar.gz

genus_extract=("{'" "'.join(self.genera_names)}")

for genus in "${{genus_extract[@]}}"; do
    echo $genus
    filename="{self.genera_table_dir}/${{genus// /_}}.tsv"
    head -n2 {self.assembly_dump_fp} | tail -n 1 | cut -c 3- > $filename
    grep "$genus" {self.assembly_dump_fp} >> "$filename"
done
"""
        call_download = f"""
python scripts/subset_and_download_genomes.py --genera-tables {' '.join(x.as_posix() for x in self.genera_table_fps)} -o {self.output_dir}

cat {self.genomes_dir}/GCA*.fna.gz > {self.output_dir}/combined_genomes.fna.gz
"""
        if self.genera_tables_exists:
            self._syscall=call_download
        else:
            self._syscall = "\n".join([call_assembly_dump, call_download])

