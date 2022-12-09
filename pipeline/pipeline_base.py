import configparser
from enum import Enum
import logging
import re
import time
from typing import List, Tuple, Union, Dict
from pathlib import Path
import os
import humanize
import pandas as pd
import graphviz

import numpy as np
from qsub_modules.base import Base as QBase

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
        config_file: Path,
        pipe_name:str = Path(__file__).stem,
        max_workers:int = 5,
        iteration_sleep:int = 120,
        rerun_downstream = False,
        testing: bool = False,
        ):

        
        
        config_path = Path(config_file)
        if not config_file.is_file():
            raise FileNotFoundError("Config file not found -> "+config_file.as_posix())
        self.project_config = configparser.ConfigParser()
        self.project_config.read(config_path)

        self.max_workers = max_workers
        self.iteration_sleep = iteration_sleep
        self.testing = testing
        self.pipe_name = pipe_name
        self.rerun_downstream = rerun_downstream

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
        def children_recursive(parent, mapping, old_children=set()):   
            new_children = mapping.get(parent, set())
            old_children.update(new_children)
            for child in new_children:
                children_recursive(child, mapping, old_children)
            return old_children


        self.log.info("Checking if any jobs are already finished")
        for label, cls in self.job_map.items():
            self.log.debug(f"Checking {label}")
            if cls.is_successful:
                self.log.debug(f"{label} has already been run - moving to completed.")
                self.jobs_holding.remove(label)
                self.jobs_complete.add(label)

        if self.rerun_downstream:
            # if any upstream dependencies have been triggered downstream should be rerun.
            jobs_to_rerun = set()
            for job_id in self.jobs_holding:
                jobs_to_rerun.update(children_recursive(job_id, self.parentage_dict))
            if jobs_to_rerun:
                self.log.warning(f"Jobs marked for rerun due to upstream changes ({len(jobs_to_rerun)})")
                self.jobs_holding.update(jobs_to_rerun)
                self.jobs_complete -= jobs_to_rerun



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
                    return {'color':'deepskyblue:dodgerblue3', 'style':'filled'}
            elif job_id in self.jobs_complete:
                return {'color':'darkolivegreen4', 'style':'filled'}
            elif job_id in self.jobs_failed:
                return {'color':'red', 'style':'filled'}
            else:
                return {'color':'grey', 'style':'filled'}


    def _add_nodes(self, dot, jobs_total) -> list:
        jobs1 = [x.split(".")[0] for x in jobs_total if len(x.split("."))==1]
        jobs2 = pd.DataFrame([x.split(".") for x in jobs_total if len(x.split("."))==2], columns=["cluster","member"])
        jobs3 = pd.DataFrame([x.split(".") for x in jobs_total if len(x.split("."))==3], columns=["cluster","subcluster","member"])
        
        clusters = set()
        # Generate nodes
        #level1 jobs
        for member_name in jobs1:
            dot.node(member_name, **self._style_node(member_name))

        #level2 jobs
        for c_name, df_c in jobs2.groupby("cluster"):
            c_label = f"cluster_{c_name}"
            clusters.add(c_name)
            with dot.subgraph(name=c_label) as c:
                c.attr(shape='box', style= 'rounded,filled', color='#d4d0cf', label=c_name)
                member_names = [".".join(x) for x in sorted(df_c.values.tolist(), key=lambda x:x[1])]
                for i in range(len(member_names)-1):
                    label = member_names[i].split(".")[-1]
                    c.node(member_names[i], label=label, **self._style_node(member_names[i]))
                    c.edge(member_names[i], member_names[i+1], style="invis", constraint="false")
                c.node(member_names[-1], label=member_names[-1].split(".")[-1], **self._style_node(member_names[-1]))

        #level3 jobs
        for (c_name, sc_name), df_c in jobs3.groupby(["cluster","subcluster"]):
            clusters.add(c_name)
            with dot.subgraph(name=f"cluster_{c_name}") as c:
                c.attr(shape='box', style= 'rounded,filled', color='#d4d0cf', label=c_name)
                clusters.add(f"{c_name}.{sc_name}")
                with c.subgraph(name=f"cluster_{c_name}.{sc_name}") as sc:
                    member_names = [".".join(x) for x in sorted(df_c.values.tolist(), key=lambda x:x[2])]
                    for i in range(len(member_names)-1):
                        label = member_names[i].split(".")[-1]
                        sc.attr(rankdir="TB")
                        sc.attr(shape='box', style= 'rounded,filled', color='#bdb5b3', label=sc_name)
                        sc.node(member_names[i], label=label, **self._style_node(member_names[i]))
                        sc.edge(member_names[i], member_names[i+1], style="invis")
                    sc.node(member_names[-1], label=member_names[-1].split(".")[-1], **self._style_node(member_names[-1]))
        return clusters

    def _generate_edge(self, child:str, parent:str, clusters:set):
        child_split = child.split(".")
        child_depth = len(child_split)
        
        parent_split = parent.split(".")
        parent_depth = len(parent_split)
        
        edge_level = min(child_depth, parent_depth)
        #print(edge_level)
        #max_shared_index = edge_level-1
        
        parent_shared_index = ".".join(parent_split[0:edge_level])
        if parent_shared_index in clusters:
            tail = "cluster_"+ parent_shared_index
        else:
            tail = parent_shared_index
        
        child_shared_index = ".".join(child_split[0:edge_level])
        
        if child.startswith("analysis"):
            head = child
        else:
            if child_shared_index in clusters:
                head = "cluster_"+ child_shared_index
            else:
                head = child_shared_index

        #if within cluster we should set constrain = False to simplify.
        constraint="true"
        if (poten_cluster := parent.rsplit(".",1)[0]) == child.rsplit(".",1)[0]:
            if poten_cluster in clusters:
                constraint = "false"
            
        return {
            'tail_name' :parent,
            'head_name':child,
            'ltail': tail,
            'lhead': head,
            'constraint':constraint 
        }
    
    def _add_edges(self, dot, dep_dict, clusters):
        all_edges = []
        for child, parents in dep_dict.items():
            for parent in parents:
                all_edges.append(self._generate_edge(child, parent, clusters))

        df_edges = pd.DataFrame(all_edges)
        df_edges.sort_values(['head_name', 'tail_name'], ascending=[True, False], inplace=True)
        df_primary_edges = df_edges.drop_duplicates(["ltail","lhead"], keep="first").reset_index(drop=True)
        
        #remove non-clusters from lhead/ltail
        df_primary_edges.loc[~df_primary_edges.loc[:,"lhead"].str.startswith("cluster_"),"lhead"] = None
        df_primary_edges.loc[~df_primary_edges.loc[:,"ltail"].str.startswith("cluster_"),"ltail"] = None
        
        df_primary_edges.to_csv("test_edges.csv", index=False)
        for edge in df_primary_edges.to_dict("records"):
            dot.edge(**edge)
        

    def generate_graphviz(self, file:Path):
        #Path("test/dep_dict.json").write_text(str(self.dependency_dict))
        #Path("test/all_jobs").write_text(str(self.jobs_total))
        dot = graphviz.Digraph(node_attr={'shape':'box', 'style': 'rounded,filled'})
        dot.attr(rankdir="LR", compound="true")
        dot.attr("node", style='filled')

        cluster_set = self._add_nodes(dot, self.jobs_total)
        self._add_edges(dot, self.dependency_dict, cluster_set)

        dot.render(file, cleanup=True)
        try:
            self.log.info(f"Dependency graph written -> {file}.pdf")
        except Exception:
            self.log.error("Failed call to DOT ensure module is loaded: [graphviz/2.40.1]")

    @property
    def log_setup(self):
        return dict(
                name = self.pipe_name,
                level = logging.getLevelName(self.project_config.get("ProjectWide","LoggingLevel")),
                log_file = (Path("logs/pipeline")/self.pipe_name).with_suffix(".log")
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

                log_file_crit = log_file.with_suffix(".critical")
                log_file_crit.unlink(missing_ok=True)
                file_handler_crit = logging.FileHandler(log_file_crit)
                file_handler_crit.setFormatter(formatter)
                file_handler_crit.setLevel(logging.WARNING)
                logger.addHandler(file_handler_crit)


                logger.debug(f"logfile at -> {log_file}")
            self._log = logger
        return self._log

    def build_dependency_dict(self):
        all_jobs = set()
        dependency_dict = {}
        parentage_dict  = {}
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
            for parent in parents:
                if parent in parentage_dict:
                    parentage_dict[parent].update(children)
                else:
                    parentage_dict[parent] = children

        self.jobs_holding = all_jobs.copy()
        self.jobs_total   = all_jobs.copy()
        self.dependency_dict = dependency_dict
        self.parentage_dict  = parentage_dict
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

                if not job_cls.is_successful:
                    self.log.warning(f"{job} failed initial success-check - retrying in (15) seconds")
                    time.sleep(15)

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
        if self.jobs_failed:
            msg = f"Failed number of jobs: ({len(self.jobs_failed)})"
            self.log.error(msg)
            raise RuntimeError(msg)