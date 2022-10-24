from functools import partial
from pathlib import Path
from re import S
from Bio import SeqIO
from Bio.Seq import Seq
import configparser
import logging
from typing import List, Literal, Sequence, Tuple, Union, Set
from dataclasses import dataclass, field


import numpy as np

from scripts.kmer_gen_catalogue import canonicalize

@dataclass
class BGCData:
    name:               str = field(default_factory=str)
    sequence:           str = field(default_factory=str, repr=False)
    origin_path:        Path = field(default_factory=Path, repr=False)

    kmers:              List[str] = field(default_factory=list, repr=False)
    kmers_unique:       Set[str] = field(default_factory=set, repr=False)

    kmers_cannon:              List[str] = field(default_factory=list, repr=False)
    kmers_unique_cannon:       Set[str] = field(default_factory=set, repr=False)

@dataclass
class BGCSuperCluster:
    name:   str = field(default_factory=str)
    cluster_type:   str = field(default_factory=str)
    members: List[BGCData] = field(default_factory=list, repr=False)
    size: int = field(default_factory=int) 

    @property
    def kmers_within(self) -> dict:
        if not hasattr(self, "_kmers_within"):
            shared_kmers = dict()
            for m in self.members:
                for k in m.kmers_unique:
                    if k in shared_kmers:
                        shared_kmers[k] += 1/self.size
                    else:
                        shared_kmers[k] = 1/self.size
            self._kmers_within = shared_kmers
        return self._kmers_within
        
    kmers_distinct: List = field(default_factory=list, repr=False)
    def describe_distinct(self):
        if self.kmers_distinct == []:
            return None
        data = np.array(list(self.kmers_distinct.values()))
        quartiles = np.percentile(data, [25, 50, 75])
        data_min, data_max = data.min(), data.max()
        outstr = f"""\
{self.name} distribution of distinct kmers.
Min:    {int(data_min*self.size)}/{self.size}
Q1:     {int(quartiles[0]*self.size)}/{self.size}
Median: {int(quartiles[1]*self.size)}/{self.size}
Q3:     {int(quartiles[2]*self.size)}/{self.size}
Max:    {int(data_max*self.size)}/{self.size}\
"""
        return outstr

class CatalogueAssembler:
    def __init__(self):
        project_config = configparser.ConfigParser()
        project_config.read("../config/project_config.ini") #TODO: can we get around this relative import?
        self.project_config = project_config
        
        self.kmerlength = project_config.getint("KmerQuantification","KmerLength")

    ## Front matter
    @property
    def log_setup(self):
        return dict(
                name = self.__class__.__name__,
                level = logging.getLevelName(self.project_config.get("ProjectWide","LoggingLevel")),
                log_file = None
            )    
    @property
    def log(self):
        if not hasattr(self, "_log"):
            setup = self.log_setup
            logger = logging.getLogger(setup["name"])
            
            #remove old handlers:
            while logger.handlers:
                logger.removeHandler(logger.handlers[0])
            
            logger.setLevel(setup["level"])
            
            F = "[%(asctime)s %(name)s:%(funcName)s]%(levelname)s: %(message)s"
            formatter = logging.Formatter(F, datefmt='%d-%b-%y %H:%M:%S')
            stream_handler = logging.StreamHandler()
            stream_handler.setFormatter(formatter)
            stream_handler.setLevel(setup["level"])
            logger.addHandler(stream_handler)
            
            if setup["log_file"]:
                Path(setup["log_file"]).parent.mkdir(parents=True, exist_ok=True)
                file_handler = logging.FileHandler(setup["log_file"])
                file_handler.setFormatter(formatter)
                file_handler.setLevel(setup["level"])
                logger.addHandler(file_handler)
                logger.debug("logfile at -> "+setup["log_file"])
            self._log = logger
        return self._log
    
    @property
    def bgcs(self) -> List[BGCData]:
        if not hasattr(self, "_bgcs"):
            msg = "No BGCs has been loaded."
            self.log.error(msg)
            raise ValueError(msg)
        return self._bgcs

    @bgcs.setter
    def bgcs(self, bgcs: Union[BGCData, List[BGCData]]):
        if isinstance(bgcs, BGCData):
            self._bgcs =  [BGCData]
        elif all(isinstance(bgc, BGCData) for bgc in bgcs):
            self._bgcs = sorted(bgcs, key=lambda bgc:bgc.name)
        else:
            raise ValueError("bgcs should be set as a list of BGCData dataclass elements")

    @property
    def bgc_names(self) -> List[str]:
        return [bgc.name for bgc in self.bgcs]

    @property
    def bgc_sequences(self) -> List[str]:
        return [bgc.seqeunce for bgc in self.bgcs]

    
    def load_bgcs_from_file(self, filepath: Union[Path, Sequence[Path]], filetype: str='fasta'):
        self.log.debug("Loading BGCS")
        if not isinstance(filepath, list):
            filepath = [filepath]
        bgcs = []
        for fp in filepath:
            with open(fp) as fh:
                for record in SeqIO.parse(fp, filetype):
                    self.log.debug(f"Loading: {record.name}")
                    bgcs.append(
                        BGCData(
                            name=record.name,
                            sequence=record.seq,
                            origin_path=fp
                        )
                    )
        self.log.info(f"Loaded ({len(bgcs)}) BGC(s)")
        self.bgcs = bgcs

    @property
    def bgc_kmers(self) -> List[List[str]]:
        if not all(bgc.kmers for bgc in self.bgcs):
            self.generate_kmers()
        return [bgc.kmers for bgc in self.bgcs]

    @property
    def bgc_kmers_unqiue(self) -> List[List[str]]:
        if not all(bgc.kmers_unique for bgc in self.bgcs):
            self.generate_kmers()
        return [bgc.kmers_unique for bgc in self.bgcs]

    @property
    def bgc_kmers_cannon(self) -> List[List[str]]:
        if not all(bgc.kmers_cannon for bgc in self.bgcs):
            self.generate_kmers()
        return [bgc.kmers_cannon for bgc in self.bgcs]

    @property
    def bgc_kmers_unique_cannon(self) -> List[List[str]]:
        if not all(bgc.kmers_unique_cannon for bgc in self.bgcs):
            self.generate_kmers()
        return [bgc.kmers_unique_cannon for bgc in self.bgcs]


    def kmerise(self, seq, k:int=21):
        if len(seq) < k:
            raise ValueError(f"k cannot be larger than lenght of sequence -> k={k} > len(seq)={len(seq)}")
        return [seq[i:i+k] for i in range(len(seq)-k+1)]

    def kmerise_random(self, sequence: str, sample_size:int=1*10**3, k:int=21):
        all_kmers = self.kmerise(sequence, k)
        return np.random.choice(all_kmers, sample_size, replace=False)

    def kmerise_spaced(self, sequence:str, k:int=21, spacing:int=10):
        i = 0
        max_i = len(sequence)
        kmers = []
        while i+k <= max_i:
            kmers.append(sequence[i:i+k])
            i += k + spacing
        return kmers

    def canonicalize(self, mer:str) -> str:
        mer1 = Seq(mer)
        mer2 = mer1.reverse_complement()
        if mer1 < mer2:
            return mer1
        return mer2

    def use_kmers_all(self):
        self._kmer_generator = partial(self.kmerise, k=self.kmerlength)

    def use_kmers_random(self, sample_size=1*10**3):
        self._kmer_generator = partial(self.kmerise_random, sample_size=sample_size, k=self.kmerlength)

    def use_kmers_spaced(self, spacing=10):
        self._kmer_generator = partial(self.kmerise_spaced, spacing=spacing, k=self.kmerlength)

    @property
    def kmer_generator(self):
        if not hasattr(self, "_kmer_generator"):
            default_strategy = partial(self.kmerise, k=self.kmerlength)
            self.log.warning(f"Using default kmer strategy. {default_strategy}")
            return default_strategy
        return self._kmer_generator

    def generate_kmers(self) -> None:
        kmer_generator = self.kmer_generator
        self.log.info(f"Generating kmers using generator: {kmer_generator}")

        for bgc in self.bgcs:
            bgc.kmers           = kmer_generator(bgc.sequence)
            bgc.kmers_cannon    = [self.canonicalize(kmer) for kmer in bgc.kmers]

            bgc.kmers_unique        = set(bgc.kmers)
            bgc.kmers_unique_cannon = set(bgc.kmers_cannon)

    @property
    def superclusters(self) -> List[BGCSuperCluster]:
        if not hasattr(self, "_superclusters"):
            msg = "Clusters not assigned - use .apply_superclusters()"
            self.log.error(msg)
            raise RuntimeError(msg)
        return self._superclusters

    def apply_superclusters(self, mapping:dict, family_type:str):
        """
        arrange bgcs in superclusters based on dict with structure:
        mapping = {
            family_name = [members, ...]
        }
        family_type = name of family grouping, ie. BiG-SCAPE or other.
        """
        superclusters = []
        for super_name, super_members in mapping.items():
            SC = BGCSuperCluster(name=super_name, type=family_type)
            for member in super_members:
                if isinstance(member, BGCData):
                    SC.members.append(member)
                    SC.size += 1
                else:
                    #Check if the member is loaded.
                    loaded = [bgc for bgc in self.bgcs if bgc.name == member]
                    if len(loaded) == 1:
                        SC.members.append(loaded[0])
                        SC.size += 1
                    else:
                        self.log.warning(f"{member} in {super_name} not found.")
            superclusters.append(SC)
        self._superclusters = superclusters

    def assign_distinct_sc_kmers(self) -> None:
        for sc in self.superclusters:
            outkmers = set()
            for sc_out in self.superclusters:
                if sc == sc_out: 
                    continue
                outkmers = outkmers.union(set(sc_out.kmers_within.keys()))
            sc.kmers_distinct = {k:v for k,v in sc.kmers_within.items() if k not in outkmers}
            self.log.debug("Frequency of distinct_kmers\n"+sc.describe_distinct())
    
    def print_catalogues(self, directory:Path):
        self.log.info(f"Writing catalogue to: {directory}")
        directory = Path(directory)

        if not directory.is_dir():
            self.log.info(f"Catalogue dir does not exist - creating...")
            directory.mkdir(parents=True)
        for sc in self.superclusters:
            catalogue_file = directory / (sc.name+".catalogue")
            catalogue_entries = [f">kmer {i}\n{mer}" for i, mer in enumerate(sc.kmers_distinct)]
            catalogue_file.write_text("\n".join(catalogue_entries))





    