from pathlib import Path
from re import S
from Bio import SeqIO
from Bio.Seq import Seq
import configparser
import logging
from typing import List, Sequence, Union, Set
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
    kmers_signature:    Set[str] = field(default_factory=set, repr=False)

    kmers_cannon:              List[str] = field(default_factory=list, repr=False)
    kmers_unique_cannon:       Set[str] = field(default_factory=set, repr=False)
    kmers_signature_cannon:    Set[str] = field(default_factory=set, repr=False)

@dataclass
class BGCSuperCluster:
    name:   str = field(default_factory=str)
    type:   str = field(default_factory=str)
    members: List[BGCData] = field(default_factory=list)

    
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

    @staticmethod
    def kmerise(seq, k:int=21):
        if len(seq) < k:
            raise ValueError(f"k cannot be larger than lenght of sequence -> k={k} > len(seq)={len(seq)}")
        return [seq[i:i+k] for i in range(len(seq)-k+1)]

    @staticmethod
    def canonicalize(mer:str) -> str:
        mer1 = Seq(mer)
        mer2 = mer1.reverse_complement()
        if mer1 < mer2:
            return mer1
        return mer2

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

    def generate_kmers(self) -> None:
        self.log.info("Generating kmers")
        for bgc in self.bgcs:
            bgc.kmers           = self.kmerise(bgc.sequence, k=self.kmerlength)
            bgc.kmers_cannon    = [self.canonicalize(kmer) for kmer in bgc.kmers]

            bgc.kmers_unique        = set(bgc.kmers)
            bgc.kmers_unique_cannon = set(bgc.kmers_cannon)

    def generate_poc_signatures(self) -> None:
        self.log.info("Generating Proof-of-Concept signature")
        if not all(bgc.kmers for bgc in self.bgcs):
            self.generate_kmers()

        for bgc in self.bgcs:
            self.log.debug(f"Generating signature for {bgc.name}")

            kmers_outgroup        = set()
            kmers_outgroup_cannon = set()
            for bgc_out in self.bgcs:
                if bgc_out == bgc:
                    continue
                kmers_outgroup = kmers_outgroup.union(bgc_out.kmers_unique)
                kmers_outgroup_cannon = kmers_outgroup_cannon.union(bgc_out.kmers_unique_cannon)

            bgc.kmers_signature         = bgc.kmers_unique.difference(kmers_outgroup)
            bgc.kmers_signature_cannon  = bgc.kmers_unique_cannon.difference(kmers_outgroup_cannon)
            self.log.debug(f"Number of signature kmers: {len(bgc.kmers_signature)}")
            self.log.debug(f"Number of signature kmers(cannonical): {len(bgc.kmers_signature_cannon)}")

    def print_catalogues(self, directory:Path):
        self.log.info(f"Writing catalogue to: {directory}")
        directory = Path(directory)
        

        if not directory.is_dir():
            self.log.info(f"Catalogue dir does not exist - creating...")
            directory.mkdir(parents=True)
        for bgc in self.bgcs:
            catalogue_file = directory / (bgc.name+".catalogue")
            catalogue_entries = [f">kmer {i}\n{mer}" for i, mer in enumerate(bgc.kmers_signature)]
            catalogue_file.write_text("\n".join(catalogue_entries))





    