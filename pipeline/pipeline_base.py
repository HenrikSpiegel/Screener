import configparser
from enum import Enum
import logging
from typing import List, Tuple, Union, Dict
from pathlib import Path
import os
from scripts.qsub_base import Base as QBase

project_config = configparser.ConfigParser()
config_path = Path(__file__).parent.parent/'config'/'project_config.ini'
project_config.read(config_path)

class QsubStatus(Enum):
    UNKNOWN = None
    RUNNING  = "R"
    HOLDING  = "H"
    QUEING   = "Q"
    COMPLETE = "C"
  

class PipelineBase:
    QsubStatus = QsubStatus
    jobs_holding    = set()
    jobs_queing     = set()
    jobs_running    = set()
    jobs_complete   = set()
    jobs_failed   = set()

    def __init__(self, 
        dependencies: List[Tuple[Union[str, List[str]], Union[str, List[str]]]],
        job_map: Dict[str, QBase],
        ):

        self.log.info("Building dependency tree")
        self.dependencies = dependencies
        self.build_dependency_dict()

        self.job_map = job_map

        self.add_logging()
        self._initial_check()
        self.log.info(f"Initial status:\n{self.status}")

    def add_logging(self):
        """
        Add pipelines log to qsub classes
        """
        self.log.debug("Grab logs")
        for cls in self.job_map.values():
            cls.add_external_log(self.log)
    
    def _initial_check(self):
        """
        Move jobs which has already been completed.
        """
        self.log.info("Checking if any jobs are already finished")
        for label, cls in self.job_map.items():
            self.log.debug(f"Checking {label}")
            if cls.successful():
                self.log.info(f"{label} has already been run - moving to completed.")
                self.jobs_holding.remove(label)
                self.jobs_complete.add(label)

    @property
    def status(self):
        status_str = f"""\
jobs_holding: ({len(self.jobs_holding)})   
{self.jobs_holding}
jobs_queing: ({len(self.jobs_queing)})
{self.jobs_queing}
jobs_running: ({len(self.jobs_running)})
{self.jobs_running}
jobs_complete: (({len(self.jobs_complete)}))
{self.jobs_complete}
jobs_failed: (({len(self.jobs_failed)}))
{self.jobs_failed}
"""
        return status_str


    @property
    def log_setup(self):
        return dict(
                name = self.__class__.__name__,
                level = logging.getLevelName(project_config.get("ProjectWide","LoggingLevel")),
                log_file = "../logs/pipeline/test.log"
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
                Path(os.path.dirname(setup["log_file"])).parent.mkdir(parents=True, exist_ok=True)
                file_handler = logging.FileHandler(setup["log_file"])
                file_handler.setFormatter(formatter)
                file_handler.setLevel(setup["level"])
                logger.addHandler(file_handler)
                logger.debug("logfile at -> "+setup["log_file"])
            self._log = logger
        return self._log

    def build_dependency_dict(self):
        all_jobs = set()
        dependency_dict = {}
        for dep in self.dependencies:
            parents = dep[0] if isinstance(dep[0], set) else {dep[0]}
            children = dep[1] if isinstance(dep[1], set) else {dep[1]}
            all_jobs.update(parents ^ children)
            for child in children:
                if child in dependency_dict:
                    dependency_dict[child].update(parents)
                else:
                    dependency_dict[child] = parents
        self.jobs_holding = all_jobs
        self.dependency_dict = dependency_dict
        self.log.info(f"Found ({len(all_jobs)}) jobs to run.")
        self.log.debug(str(all_jobs))

    def check_running_jobs(self):
        self.log.warning("No running jobs")
        for job in self.jobs_running:
            if not hasattr(job, "job_id") or job.job_id=="NoID":
                self.log.error("Job without ID -> must have failed somehow. Moving to failed.")
            if job.qstat_status != self.QsubStatus.COMPLETE:
                self.log.debug('{job} not finished')

    def job_has_upstream_dependencies(self, job):
        self.log.debug(f"Checking upstream dependencies for {job}")
        upstream_jobs = self.dependency_dict.get(job, {})
        has_upstream = False

        if not upstream_jobs:
            return has_upstream
        
        for upstream_job in upstream_jobs:
            if upstream_job in self.jobs_complete:
                self.log.debug(f"--Complete upsteam: {upstream_job}")
            elif upstream_job in self.jobs_holding:
                self.log.debug(f"--Holding upsteam: {upstream_job}")
                has_upstream = True
            elif upstream_job in self.jobs_queing:
                self.log.debug(f"--Queing upstream: {upstream_job}")
                has_upstream = True
            elif upstream_job in self.jobs_running:
                self.log.debug(f"--Running upstream: {upstream_job}")
                has_upstream = True
            
        return has_upstream

    def move_to_que(self, ids:set):
        if isinstance(ids, str): 
            ids = {ids}
        self.jobs_queing.update(ids)
        self.jobs_holding -= ids

    def move_to_running(self, ids:set):
        if isinstance(ids, str): 
            ids = {ids}
        self.jobs_running.update(ids)
        self.jobs_queing -= ids

    def move_to_complete(self, ids:set):
        if isinstance(ids, str): 
            ids = {ids}
        self.jobs_complete.update(ids)
        self.jobs_running -= ids
    
    def find_ready_jobs(self):
        ready_jobs = set()
        for job_id in self.jobs_holding:
            # Check if has dependencies:
            if self.job_has_upstream_dependencies(job_id):
                continue
            self.log.debug(f"{job_id}: No upstream dependencies.")
            ready_jobs.add(job_id)
        self.log.info(f"Added ({len(ready_jobs)}) jobs to que")
        self.log.debug(f"Adding jobs to que -> {ready_jobs}")
        self.move_to_que(ready_jobs)



    
    
    