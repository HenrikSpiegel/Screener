import logging
import os, sys
import pathlib
from scripts.functions import submit2

class Base:
    working_dir = "/home/projects/dtu_00009/people/henspi/git/Screener"

    #Configs: should be moved to a real config file:
    @property
    def working_dir(self):
        return "/home/projects/dtu_00009/people/henspi/git/Screener"

    def as_abspath(self, path):
        #simple wrapper to add the wdir to a relative path.
        return os.path.join(self.working_dir, path)
    
    @property
    def log_setup(self):
        return dict(
                name = self.__class__.__name__,
                level = logging.DEBUG,
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
                pathlib.Path(os.path.dirname(setup["log_file"])).mkdir(parents=True, exist_ok=True)
                file_handler = logging.FileHandler(setup["log_file"])
                file_handler.setFormatter(formatter)
                file_handler.setLevel(setup["level"])
                logger.addHandler(file_handler)
                logger.debug("logfile at -> "+setup["log_file"])
            self._log = logger
        return self._log
    
    def add_external_log(self, log: logging.Logger):
        self.log.addHandler(log.handlers[1])
        # self.log.debug()
        # self.
        # for handler in log.handlers:
        #     if isinstance(handler, logging.StreamHandler):
        #         continue
        #     self.log.debug(f"Adding -> {handler}")
        #     self._log.addHandler(handler)

  
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

    def set_qsub_args(self, jobtag: str=None, **kwargs):
        jobname = self.__class__.__name__
        if jobtag:
            jobname = "_".join([jobname, jobtag])

        if not hasattr(self, "qsub_requirements"):
            self.log.warning("Uses qsub requirements from qsub_base - set 'qsub_requirements' property (dict)")
            qsub_args = self.default_qsub_requirement
        else:
            qsub_args = self.qsub_requirements

        qsub_args.update(dict(
            directory = self.working_dir,
            group="dtu_00009",
            jobname=jobname,
            output = os.path.join("logs", "qsub", jobname+ "_stdout"),
            error = os.path.join("logs","qsub", jobname+ "_stderr")
        ))
        qsub_args.update(kwargs)
        self._qsub_args = qsub_args
    
    @property
    def syscall(self):
        if not hasattr(self, "_syscall"):
            self.log.warning("running default syscall")
            self.generate_syscall()
        return self._syscall

    def generate_syscall(self, **kwargs):
        raise NotImplementedError
        
    def preflight(self, **kwargs) -> None:
        raise NotImplementedError
        
    def output_exist(self, **kwargs) -> bool:
        raise NotImplementedError
    
    def is_succes(self) -> bool:
        raise NotImplementedError

    @staticmethod
    def gen_prefix(reads_gb:int) -> str:
        return str(reads_gb).replace(".","_")+"GB"

    def add_to_que(self, test=False) -> None:
        self.log.debug(f"command:\n{self.syscall}")
        id = submit2(
            command = self.syscall,
            **self.qsub_args,
            test=test
        )
        self.log.info(f"Added qsub job: {self.qsub_args['jobname']} / {id}")
        self.job_id = id

