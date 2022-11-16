from qsub_modules.base import Base
import os
import pathlib
from pathlib import Path
import glob
import logging
from typing import List, Union

class QuantifierKmer(Base):
    def __init__(self, read_files:Union[Path, List[Path]], fp_catalogue:Union[Path, List[Path]], output_dir: Path, kmer_size = 21, log: logging.Logger=None) -> None:
        if log:
            self.add_external_log(log)
        
        if isinstance(read_files, str):
            read_files = [read_files]
        self.read_files      = [Path(x) for x in read_files]
        self.fp_catalogue = Path(fp_catalogue)

        output_dir = Path(output_dir)
        self.output_dir = output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.success_file = output_dir / "success"

        self.kmer_size = kmer_size

        self.fp_countdir = output_dir/"counts"
        self.fp_countdir.mkdir(parents=True, exist_ok=True)
        self.fp_readmer = output_dir / f"read_{kmer_size}mers_.jf"

    qsub_requirements = dict(
        modules = "tools anaconda3/2021.05 jellyfish/2.3.0",
        runtime = 60,
        cores = 20,
        ram=60,
        )
    
    # @property
    # def log_setup(self):
    #     setup = super().log_setup
    #     setup.update({"level":logging.INFO})
    #     return setup

    def preflight(self, check_input=False) -> None:
        fp_catalogue = self.fp_catalogue 
        if isinstance(fp_catalogue, list):
            self.catalogues  = [Path(x) for x in fp_catalogue]
        else:
            fp_catalogue = Path(fp_catalogue)
            if fp_catalogue.is_dir():
                catalogues = list(fp_catalogue.glob("*.catalogue"))
                self.log.debug(f"Catalogue size: ({len(catalogues)})")
                self.catalogues  = catalogues
            elif fp_catalogue.is_file():
                self.catalogues  = catalogues
            else:
                self.log.error("Failed parsing catalogue")
                raise RuntimeError("Failed parsing catalogue")

        self.fp_catalogue_index = self.output_dir/"catalogue.index"
        bgc_name = [x.stem for x in self.catalogues]

        
        lines = [name+","+path.as_posix() for name,path in zip(bgc_name, self.catalogues)]
        self.fp_catalogue_index.write_text(("\n".join(lines))+"\n") #For some reason the last line is ignored unless there is a newline.


        if check_input:
            catalogues_exist = [os.path.isfile(x) for x in self.catalogues]
            if not all(catalogues_exist):
                err_msg = f"Catalogue files not found -> {[file for file, exists in zip(self.catalogues, catalogues_exist) if not exists]}"
                self.log.error(err_msg)
                raise IOError(err_msg)

            
            if not all(os.path.isfile(file) for file in self.read_files):
                self.log.error("Reads file not found -> "+self.read_files)
                raise IOError("Reads file not found -> "+self.read_files)
            self.log.debug("Found all input files")

    def generate_syscall(self) -> None:
        # sets mem / cpu based on default qsub args.

        # if not self.reference.endswith(".gbk"):
        syscall=f"""
set -e
# Create kmer database
echo "Creating DB"
zcat {" ".join(x.as_posix() for x in self.read_files)} | jellyfish count -m {self.kmer_size} -s 5G -t {self.qsub_args['cores']-1} -C -o {self.fp_readmer} /dev/fd/0

# Get counts for catalogue(s)
in="${{1:-{self.fp_catalogue_index}}}"
[ ! -f "$in" ] && {{ echo "$0 - File $in not found."; exit 1; }}
while IFS= read -r file
do
    name="$(cut -d',' -f1 <<<$file)"
    path="$(cut -d',' -f2 <<<$file)"
	## avoid commented filename ##
	#[[ $file = \#* ]] && continue
    call="jellyfish query {self.fp_readmer} -s $path"
    echo "Counting for: $name"
    echo "\t($path)"
	$call > {self.fp_countdir}/$name.counted
done < "${{in}}"

# Summarise:
average_readlength=$(python scripts/average_readlength.py -f {" ".join(x.as_posix() for x in self.read_files)})
echo "Average readlength: $average_readlength"

python scripts/kmer_summarise.py --directory {self.fp_countdir} -o {self.output_dir / "kmer_summation.tsv"} -l $average_readlength

"""
        self._syscall = syscall
        return

    # def successful(self):
    #     success_file = self.output_dir / "success"
    #     self.log.debug("Run was succesful")
    #     return os.path.isfile(success_file)

    # @staticmethod
    # def is_success(output_dir) -> str:
    #     return os.path.isfile(os.path.join(output_dir, "kmer_summation.tsv"))




if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--reads", required=True, help="file(s) containing reads to be mapped")
    parser.add_argument("--catalogue", required=True, help="Path to catalogue file or dir with multiple catalogue files (.fa) - must exist on init.")
    parser.add_argument("-o", required=True, help="output dir")
    parser.add_argument("-k", default=21, help="Kmer length [21]")
    args = parser.parse_args()

    api = QuantifierKmer(read_files=args.reads, fp_catalogue=args.catalogue, output_dir=args.o, kmer_size = args.k)
    api.preflight(check_input=True)
    api.set_qsub_args(jobtag="test")
    api.generate_syscall() #not needed as we run the default
    api.add_to_que(test=False)