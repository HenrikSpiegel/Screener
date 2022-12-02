import subprocess
import logging
import os, sys
from pathlib import Path
import time
from qsub_modules.functions import submit2
import configparser
from enum import Enum

class QsubStatus(Enum):
    UNKNOWN = None
    RUNNING  = "R"
    HOLDING  = "H"
    QUEING   = "Q"
    COMPLETE = "C"


class Base:
    #working_dir = "/home/projects/dtu_00009/people/henspi/git/Screener"
    QsubStatus = QsubStatus
    #Configs: should be moved to a real config file:
    #@property
    #def working_dir(self):
    #    return "/home/projects/dtu_00009/people/henspi/git/Screener"
    loglvl = "DEBUG"

    def as_abspath(self, path):
        #simple wrapper to add the wdir to a relative path.
        return os.path.join(self.working_dir, path)
    
    def set_config(self, config_file:Path):
        fp = Path(config_file)
        if not fp.is_file():
            raise FileNotFoundError("Config_file doesn't exist" + fp.as_posix())
        self.config = configparser.ConfigParser()
        self.config_file = fp

    @property
    def log_setup(self):     
        return dict(
                name = self.__class__.__name__,
                level = logging.getLevelName(self.loglvl),
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
                Path(os.path.dirname(setup["log_file"])).mkdir(parents=True, exist_ok=True)

                file_handler = logging.FileHandler(setup["log_file"])
                file_handler.setFormatter(formatter)
                file_handler.setLevel(setup["level"])
                logger.addHandler(file_handler)
                logger.debug("logfile at -> "+setup["log_file"])
            self._log = logger
        return self._log
    
    def add_external_log(self, log: logging.Logger):
        self._log = log.getChild(self.log_setup["name"])

  
    @property
    def default_qsub_requirement(self) -> dict:
        args = dict(
            modules = "tools anaconda3/2021.11 ",
            runtime = 30,
            cores = 5,
            ram=20,
        )
        return args
    
    @property
    def qsub_args(self):
        if not hasattr(self, "_qsub_args"):
            self.log.warning("Using default args")
            self.set_qsub_args()
        return self._qsub_args

    def set_qsub_args(self, jobtag: str=None, jobname=None, **kwargs):
        jobname = jobname or self.__class__.__name__
        if jobtag:
            jobname = "_".join([jobname, jobtag])

        if not hasattr(self, "qsub_requirements"):
            self.log.warning("Uses qsub requirements from qsub_base - set 'qsub_requirements' property (dict)")
            qsub_args = self.default_qsub_requirement
        else:
            qsub_args = self.qsub_requirements

        qsub_args.update(dict(
            #directory = self.working_dir,
            group="dtu_00009",
            jobname=jobname,
            output = Path("logs")/"qsub"/ (jobname+ "_stdout"),
            error = Path("logs")/"qsub"/ (jobname+ "_stderr")
        ))
        qsub_args.update(kwargs)
        self._qsub_args = qsub_args
    
    @property
    def syscall(self):
        if not hasattr(self, "_syscall"):
            self.generate_syscall()
        return self._syscall

    @property
    def wrapped_syscall(self):
        call_raw = self.syscall
        call_wrapped = f"""\
main () {{
{call_raw}
}}
time main

touch {self.success_file}
sleep 10
"""
        return call_wrapped

    @property
    def qstat_dict(self) -> dict:
        if not hasattr(self, "job_id") or not self.job_id:
            self.log.error("No job id - first add to que.")
            return {}

        qstat_return = subprocess.run(["qstat", "-f", str(self.job_id)], capture_output=True)
        qstat = qstat_return.stdout.decode()
        qstat_dict = {}
        for line in qstat.split("\n"):
            line=line.strip()
            if line.startswith("exit_status"):
                qstat_dict['exit_status'] = line.split("=")[1].strip()
            elif line.startswith("job_state"):
                qstat_dict["job_state"]   = line.split("=")[1].strip()    
        return qstat_dict
    @property
    def qstat_status(self):
        status = self.qstat_dict.get("job_state", None)
        return self.QsubStatus(status)
    @property
    def qstat_exitcode(self):
        code = self.qstat_dict.get("exit_status", None)
        return code

    def generate_syscall(self, **kwargs):
        raise NotImplementedError
        
    def preflight(self, **kwargs) -> None:
        raise NotImplementedError
        
    @property
    def is_running(self):
        return self.qstat_status == QsubStatus.RUNNING

    @property
    def is_queing(self):
        return self.qstat_status == QsubStatus.QUEING

    @property
    def is_complete(self):
        return self.qstat_status == QsubStatus.COMPLETE 
    
    @property
    def is_successful(self):
        if hasattr(self, "success_file"):
            self.log.debug(self.success_file)
            return self.success_file.exists()
            # if not self.success_file.exists():
            #     self.log.warning("Initial success-check failed waiting 15 seconds")
            #     time.sleep(15)
            #     self.success_file.exists()
            # return True
        else:
            self.log.warning(f"{self.__class__.__name__} doesn't have a success_file set.")
            return False
  
    @staticmethod
    def gen_prefix(reads_gb:float) -> str:
        # returns stringified float with GB suffix
        # 1 -> 1_0GB
        # 0.1 -> 0_1GB
        return str(float(reads_gb)).replace(".","_")+"GB"

    def add_to_que(self, test=False) -> None:
        self.qsub_args["output"].parent.mkdir(parents=True, exist_ok=True)
        self.qsub_args["output"].unlink(missing_ok=True)
        self.qsub_args["error"].unlink(missing_ok=True)

        self.success_file.unlink(missing_ok=True)

        self.log.debug(f"command:\n{self.syscall}")
        id = submit2(
            command = self.wrapped_syscall,
            **self.qsub_args,
            test=test
        )
        time.sleep(0.5)
        self.log.info(f"Added qsub job: {self.qsub_args['jobname']} / {id}")
        self.job_id = id

    def wait_for_finish(self, wait_time:int=30):
        if not hasattr(self, "job_id") or self.job_id=="NoID":
            raise RuntimeError("No job_id -> nothing to wait for")
        while self.qstat_status != self.QsubStatus.COMPLETE:
            self.log.debug(f"Waiting for ({self.__class__.__name__}) to complete: ({self.qstat_status})")
            time.sleep(wait_time)
        self.log.info(f"Successful:{self.successful()} - Exit-code: {self.qstat_exitcode}")


