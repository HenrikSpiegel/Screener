
import shutil
import time
from qsub_modules.base import Base, QsubStatus
from Bio import SeqIO
from typing import Union, List
from pathlib import Path

class MCLClustering(Base):
    def __init__(self, blast_file: Path, output_dir: Path, log=None, loglvl = "DEBUG"):

        if log:
            self.add_external_log(log)
        self.loglvl = loglvl

        self.blast_file = Path(blast_file)

        self.output_dir = Path(output_dir)

        self.success_file = self.output_dir/'.success'
        self.abc_columns = "1,2,12" #qaccver, saccver and evalue



    qsub_requirements = dict(
        modules = "tools mcl/14-137", # anaconda3/2021.05 
        runtime = 360,
        cores = 20,
        ram=180,
        )

    def preflight(self, check_input=True):
        self.output_dir.mkdir(parents=True, exist_ok=True)

        if check_input and (not self.blast_file.exists()):
            raise FileNotFoundError(self.blast_file)


    def generate_syscall(self):
        
        self._syscall = f"""\
# convert blast to abc format (seq1, seq2, score(hereEvalue)) http://micans.org/mcl/man/clmprotocols.html#blast
# Also skipping the header line
cut -f {self.abc_columns} {self.blast_file} | tail -n +2  > {self.output_dir / "blast_result.abc"}

#load using mcxload to create networkfile (.mci) and dictionary file (.tab)
mcxload -abc {self.output_dir / "blast_result.abc"} --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o {self.output_dir / "blast_result.mci"} -write-tab {self.output_dir / "blast_result.tab"}

# Run the clustering
mcl {self.output_dir / "blast_result.mci"} -I 1.4  -use-tab {self.output_dir / "blast_result.tab"} -odir {self.output_dir} -te {self.qsub_requirements["cores"]-1}
mcl {self.output_dir / "blast_result.mci"} -I 2  -use-tab {self.output_dir / "blast_result.tab"} -odir {self.output_dir} -te {self.qsub_requirements["cores"]-1}
mcl {self.output_dir / "blast_result.mci"} -I 4  -use-tab {self.output_dir / "blast_result.tab"} -odir {self.output_dir} -te {self.qsub_requirements["cores"]-1}
mcl {self.output_dir / "blast_result.mci"} -I 6  -use-tab {self.output_dir / "blast_result.tab"} -odir {self.output_dir} -te {self.qsub_requirements["cores"]-1}

# collect-to-jsons.
python scripts/mcl_conv_json.py --mcl-files {self.output_dir}/*mci.I*
"""

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--blast", required=True, type=Path, help="blastoutput")
    parser.add_argument("-o", required=True, type=Path, help="output dir")
    args = parser.parse_args()

    pb = MCLClustering(
        blast_file= args.blast, 
        output_dir=args.o)
    pb.preflight()
    pb.add_to_que(test=True)


    
