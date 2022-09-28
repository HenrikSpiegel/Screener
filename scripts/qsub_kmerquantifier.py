from scripts.qsub_base import Base
import os
import pathlib
import glob
import logging
from typing import List, Union

class QuantifierKmer(Base):
    def __init__(self, read_files:Union[str, List[str]], fp_catalogue:Union[str, List[str]], output_dir: str, kmer_size = 12, log: logging.Logger=None) -> None:
        if log:
            self.add_external_log(log)
        
        if isinstance(read_files, str):
            read_files = [read_files]
        self.read_files      = read_files

        # Create a mapping between catalogue and their fps
        if isinstance(fp_catalogue, list):
            self.catalogues  = fp_catalogue
        elif os.path.isdir(fp_catalogue):
            catalogues = glob.glob(os.path.join(fp_catalogue, "*.catalogue*.fa"))
            self.log.info(f"Catalogue size: ({len(catalogues)})")
            self.catalogues  = catalogues
        else:
            self.catalogues  = [fp_catalogue]

        self.fp_catalogue_index = os.path.join(output_dir, "catalogue.index")
        bgc_name = [x.rsplit("/",1)[1].split(".catalogue")[0] for x in self.catalogues]

        pathlib.Path(output_dir).mkdir(parents=True, exist_ok=True)
        with open(self.fp_catalogue_index, "w") as fh:
            lines = [name+" "+path for name,path in zip(bgc_name, self.catalogues)]
            fh.write("\n".join(lines)+"\n")


        self.kmer_size = kmer_size
        self.output_dir = output_dir
        self.fp_countdir = os.path.join(output_dir, "counts")
        pathlib.Path(self.fp_countdir).mkdir(parents=True, exist_ok=True)
        self.fp_readmer = os.path.join(output_dir, f"read_{kmer_size}mers_.jf")

    qsub_requirements = dict(
        modules = "tools anaconda3/2021.05 jellyfish/2.3.0",
        runtime = 30,
        cores = 19,
        ram=80,
        )
    
    # @property
    # def log_setup(self):
    #     setup = super().log_setup
    #     setup.update({"level":logging.INFO})
    #     return setup

    def preflight(self, check_input=False) -> None:
        if check_input:
            catalogues_exist = [os.path.isfile(x) for x in self.catalogues]
            if not all(catalogues_exist):
                err_msg = f"Catalogue files not found -> {[file for file, exists in zip(self.catalogues, catalogues_exist) if not exists]}"
                self.log.error(err_msg)
                raise IOError(err_msg)

            
            if not all(os.path.isfile(file) for file in self.read_files):
                self.log.error("Reads file not found -> "+self.read_files)
                raise IOError("Reads file not found -> "+self.read_files)
            self.log.info("Found all input files")

        pathlib.Path(os.path.dirname(self.output_dir)).mkdir(parents=True, exist_ok=True)

    def generate_syscall(self) -> None:
        # sets mem / cpu based on default qsub args.

        # if not self.reference.endswith(".gbk"):
        syscall=f"""
set -e
# Create kmer database
echo "Creating DB"
zcat {" ".join(self.read_files)} | jellyfish count -m {self.kmer_size} -s 5G -t {self.qsub_args['cores']-1} -C -o {self.fp_readmer} /dev/fd/0

# Get counts for catalogue(s)
in="${{1:-{self.fp_catalogue_index}}}"
[ ! -f "$in" ] && {{ echo "$0 - File $in not found."; exit 1; }}
while IFS= read -r file
do
    name="$(cut -d' ' -f1 <<<$file)"
    path="$(cut -d' ' -f2 <<<$file)"
	## avoid commented filename ##
	#[[ $file = \#* ]] && continue
    call="jellyfish query {self.fp_readmer} -s $path"
    echo "Counting for: $name"
	$call > {self.fp_countdir}/$name.counted
done < "${{in}}"

# Summarise:
python -m scripts.kmer_summarise --directory {self.fp_countdir} -o {os.path.join(self.output_dir, "kmer_summation.tsv")}

# compress countfiles
gzip {self.fp_countdir}
    """
        self._syscall = syscall
        return


    @staticmethod
    def is_success(output_dir) -> str:
        return os.path.isfile(os.path.join(output_dir, "kmer_summation.tsv"))

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--reads", required=True, help="file(s) containing reads to be mapped")
    parser.add_argument("--catalogue", required=True, help="Path to catalogue file or dir with multiple catalogue files (.fa) - must exist on init.")
    parser.add_argument("-o", required=True, help="output dir")
    parser.add_argument("-k", default=21, help="Kmer length [21]")
    args = parser.parse_args()

    api = QuantifierKmer(reads=args.reads, fp_catalogue=args.catalogue, output_dir=args.o, kmer_size = args.k)
    api.preflight(check_input=True)
    api.set_qsub_args(jobtag="test")
    api.generate_syscall() #not needed as we run the default
    #print(api.syscall)
    api.add_to_que(test=False)