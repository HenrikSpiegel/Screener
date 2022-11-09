import configparser
from enum import Enum
import logging
import re
import time
from typing import List, Tuple, Union, Dict
from pathlib import Path
import os
import humanize

import graphviz

import numpy as np
from qsub_modules.base import Base as QBase

project_config = configparser.ConfigParser()
config_path = Path('config/project_config.ini')
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
    jobs_failed     = set()

    progress_updated = False

    def __init__(self, 
        dependencies: List[Tuple[Union[str, List[str]], Union[str, List[str]]]],
        job_map: Dict[str, QBase],
        pipe_name:str = Path(__file__).stem,
        max_workers:int = 5,
        iteration_sleep:int = 120,
        testing: bool = False
        ):

        self.max_workers = max_workers
        self.iteration_sleep = iteration_sleep
        self.testing = testing
        self.pipe_name = pipe_name

        self.log.info("Building dependency tree")
        self.dependencies = dependencies
        self.build_dependency_dict()

        self.job_map = job_map
        self.add_logging_from_jobs()

        self._initial_check()
        
        self.log.info(f"Initial status:\n{self.status}")
        self.generate_graphviz(self.log_setup["log_file"].as_posix()+'.gv')

    def add_logging_from_jobs(self):
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
            if cls.is_successful:
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
jobs_complete: ({len(self.jobs_complete)})
{self.jobs_complete}
jobs_failed: ({len(self.jobs_failed)})
{self.jobs_failed}
"""
        return status_str

    def _style_node(self, job_id):
            if job_id in self.jobs_holding:
                return {'color':'cornsilk', 'style':'filled'}
            elif job_id in self.jobs_queing:
                return {'color':'deepskyblue', 'style':'filled'}
            elif job_id in self.jobs_running:
                if self.job_map[job_id].is_running:
                    return {'color':'dodgerblue3', 'style':'filled'}
                else:
                    return {'color':'dodgerblue3', 'style':'striped'}
            elif job_id in self.jobs_complete:
                return {'color':'darkolivegreen4', 'style':'filled'}
            elif job_id in self.jobs_failed:
                return {'color':'red', 'style':'filled'}
            else:
                return {'color':'grey', 'style':'filled'}

    def generate_graphviz(self, file:Path):    
        
            
        sample_regex = re.compile(r"\.\d$")
        dict_belong = dict()
        dict_clusters = dict()
        for job in self.jobs_total:
            cleaned_name = sample_regex.sub("", job)
            
            if cleaned_name != job:
                clust_name = "cluster_"+cleaned_name
                dict_belong[job] = clust_name
                if clust_name in dict_clusters:
                    dict_clusters[clust_name].append(job)
                else:
                    dict_clusters[clust_name] = [job]
            else:
                dict_belong[job] = cleaned_name

        dot = graphviz.Digraph(comment=self.pipe_name, node_attr={'shape':'box', 'style': 'rounded,filled'})
        dot.attr(rankdir="LR")
        dot.attr("node", style='filled')
        clustered_nodes = set()
        for name, members in dict_clusters.items():
            clustered_nodes.update(members)
            with dot.subgraph(name="cluster_"+name) as c:
                c.attr(rankdir="TB")
                c.attr(shape='box', style= 'rounded,filled', color='#82d2f0', label=name)
                for i in range(len(members)-1):
                    sort_mem = sorted(members)
                    c.node(members[i], **self._style_node(members[i]))
                    c.edge(sort_mem[i], sort_mem[i+1], style="invis")
                c.node(members[-1], **self._style_node(members[-1]))

        for job_id in self.jobs_total - clustered_nodes:
            dot.node(job_id, **self._style_node(job_id))

        for child, parents in self.dependency_dict.items():
            if dict_belong[child].startswith('cluster_'):
                if not child == min(dict_clusters[dict_belong[child]]):
                    continue
                head = dict_belong[child]
            else:
                head = None
            for parent in parents:
                if dict_belong[parent].startswith('cluster_'):
                    if not parent == max(dict_clusters[dict_belong[parent]]):
                        continue
                    tail = dict_belong[parent]
                else:
                    tail = None
                dot.edge(parent, child, ltail=tail, lhead=head)

        dot.render(file, cleanup=True)
        self.log.info(f"Dependency graph written -> {file}.pdf")

    @property
    def log_setup(self):
        return dict(
                name = self.pipe_name,
                level = logging.getLevelName(project_config.get("ProjectWide","LoggingLevel")),
                log_file = Path("logs/pipeline")/self.pipe_name
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
                log_file = Path(setup["log_file"])
                log_file.parent.mkdir(parents=True, exist_ok=True)
                log_file.unlink(missing_ok=True)

                file_handler = logging.FileHandler(log_file)
                file_handler.setFormatter(formatter)
                file_handler.setLevel(setup["level"])
                logger.addHandler(file_handler)
                logger.debug(f"logfile at -> {log_file}")
            self._log = logger
        return self._log

    def build_dependency_dict(self):
        all_jobs = set()
        dependency_dict = {}
        for dep in self.dependencies:
            if len(dep) == 1: # use this approach ("jobname", ) to add jobs without prior dependencies
                all_jobs.add(dep[0])
                continue
            parents = dep[0] if isinstance(dep[0], set) else {dep[0]}
            children = dep[1] if isinstance(dep[1], set) else {dep[1]}
            all_jobs.update(parents ^ children)
            for child in children:
                if child in dependency_dict:
                    dependency_dict[child].update(parents)
                else:
                    dependency_dict[child] = parents
        self.jobs_holding = all_jobs.copy()
        self.jobs_total   = all_jobs.copy()
        self.dependency_dict = dependency_dict
        self.log.info(f"Found ({len(all_jobs)}) jobs to run.")
        self.log.debug(str(all_jobs))


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

    def job_has_upstream_failures(self, job):
        upstream_jobs = self.dependency_dict.get(job, {})
        if not upstream_jobs:
            return False
        for upstream_job in upstream_jobs:
            if upstream_job in self.jobs_failed:
                return True
        return False

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
        jobs_holding_iter = self.jobs_holding.copy()
        for job_id in jobs_holding_iter:
            # Check if has dependencies:
            if self.job_has_upstream_failures(job_id):
                self.progress_updated = True
                self.log.warning(f"{job_id} has upstream failures -> moving to jobs_failed")
                self.jobs_holding.remove(job_id)
                self.jobs_failed.add(job_id)
                self.progress_updated = True
                continue

            if self.job_has_upstream_dependencies(job_id):
                continue

            self.log.debug(f"{job_id}: No upstream dependencies.")
            self.log.debug(f"Adding jobs to que -> {job_id}")
            self.jobs_holding.remove(job_id)
            self.jobs_queing.add(job_id)
            self.progress_updated = True

    def check_running_jobs(self):
        if not self.jobs_running:
            self.log.warning("No jobs are running")
            return

        job_running_iterator = self.jobs_running.copy()
        for job in job_running_iterator:
            job_cls = self.job_map[job]
            self.log.debug(f"cls")
            if not hasattr(job_cls, "job_id") or job_cls.job_id=="NoID":
                self.progress_updated = True
                self.log.error("Job without ID -> must have failed somehow. Moving to failed.")
                self.jobs_failed.add(job)
                self.jobs_running.remove(job)
                continue

            if job_cls.is_running or job_cls.is_queing:
                continue

            if job_cls.is_complete:
                self.progress_updated = True
                time.sleep(2) #It appears we sometimes miss the successfile - perhaps a slight desync issue.
                if job_cls.is_successful:
                    self.log.info(f"{job} has succesfully completed")
                    self.jobs_running.remove(job)
                    self.jobs_complete.add(job)
                else:
                    self.log.error(f"{job} has failed")
                    self.jobs_running.remove(job)
                    self.jobs_failed.add(job)
                
                continue
            self.log.error(f"""\
{job} -> {job_cls}:
status: {job_cls.qstat_status}
is_running: {job_cls.is_running}
is_complete: {job_cls.is_complete}
is_successful: {job_cls.is_successful}
""")
            raise RuntimeError(f"Reached fail state during check - review code -> {job}")
    
    @property
    def jobs_to_finalize(self):
        return self.jobs_holding ^ self.jobs_queing ^ self.jobs_running

    @property
    def jobs_finalized(self):
        return self.jobs_complete ^ self.jobs_failed

    @property
    def progress(self):
        return f"{len(self.jobs_finalized)}/{len(self.jobs_to_finalize ^ self.jobs_finalized)}"
        
    @property
    def pipeline_is_finished(self):
        return self.jobs_to_finalize == set()

    @property
    def pipeline_is_stuck(self):
        """
        if nothing is running - and nothing can run we are stuck.
        """
        if self.pipeline_is_finished or (self.jobs_running ^ self.jobs_queing):
            return False
        return True

    def runjob(self, job_id:str):
        self.log.info(f"Running {job_id}")
        cls = self.job_map[job_id]
        cls.preflight(check_input = True)
        cls.set_qsub_args(jobname=job_id)
        cls.generate_syscall()
        cls.add_to_que(test = self.testing)
        try:
            int(cls.job_id)
        except Exception as err:
            raise RuntimeError(f"Qsub jobid return: {cls.job_id}")
        self.progress_updated = True

    def run_pipeline(self):
        pipeline_start = time.time()

        self.log.info("Started pipeline run.")
        self.log.info(f"Intial progress {self.progress} jobs")

        status_time_delta = 2*60
        last_status_time = time.time()

        while not self.pipeline_is_finished:
            self.check_running_jobs()
            self.find_ready_jobs()

            if self.pipeline_is_stuck:
                self.log.error(f"Stuck pipeline. Both que and running is empty (TERMINATING).\n{self.status}")
                raise RuntimeError("Pipeline is stuck")

            while len(self.jobs_running) < self.max_workers:
                self.log.debug("Unused workers")
                if not self.jobs_queing:
                    self.log.debug("No ready jobs to run")
                    break
                
                job_to_start = self.jobs_queing.pop()
                try:
                    self.runjob(job_to_start)
                except Exception as err:
                    self.log.exception(f"Failed starting {job_to_start}.")
                    self.jobs_failed.add(job_to_start)
                else:
                    self.jobs_running.add(job_to_start)

            if self.progress_updated:
                self.generate_graphviz(self.log_setup["log_file"].as_posix()+'.gv')
                self.progress_updated = False
                
            if self.pipeline_is_finished:
                break

            if (time.time() - last_status_time) > status_time_delta:
                last_status_time = time.time()
                self.log.info(f"Progress: {self.progress}\n{self.status} jobs")
                self.generate_graphviz(self.log_setup["log_file"].as_posix()+'.gv')
                
            time.sleep(self.iteration_sleep)
        
        
        self.log.info(f"Final status:\n{self.status}")
        self.generate_graphviz(self.log_setup["log_file"].as_posix()+'.gv')
        runtime = time.time() - pipeline_start
        self.log.info(f"Pipeline terminated - runtime: {humanize.naturaldelta(runtime)}")