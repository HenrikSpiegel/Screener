import argparse
from qsub_modules.base import Base
from pathlib import Path
import hashlib
import configparser
import logging


project_config = configparser.ConfigParser()
project_config.read("config/project_config.ini")


class AddToQue(Base):
    def __init__(self, command:str, success_file:str, name:str=None):
        self.command = command
        self.name = name or "AddToQue"
        self.jobtag = name

        self.success_file = success_file

    default_args = dict(
        modules = "tools anaconda3/2021.11",
        runtime = 60,
        cores = 10,
        ram=40,
        group="dtu_00009",
    )

    @property
    def log_setup(self):
        return dict(
                name = self.name,
                level = logging.getLevelName(project_config.get("ProjectWide","LoggingLevel")),
                log_file = None
            )    

    def generate_syscall(self):
        call = self.command + "\n" + f"touch {self.success_file}"
        self._syscall = call

    def preflight(self, **kwargs) -> None:
        self.log.debug(f"preflight not defined for {self.__class__.__name__}")
        pass


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--command", required=True, help="command(s) job be executed by qsub script.")
    parser.add_argument("--success", required=True, type=Path, help="location used to store successfile.")
    parser.add_argument("--name", default="addQue_def", help="jobname")
    parser.add_argument("--test", action="store_true", help="Prints the qsub job-script to stdout and doesn't add_to_que")
    parser.add_argument("--options", nargs="+", help="options in kw:arg kwargs passed to submit2 usage: -o runtime:360 cores:30 modules:'samtools/1.14 blastn/XXX'")
    args = parser.parse_args()

    api = AddToQue(
        command = args.command,
        success_file=args.success,
        test = args.test,
        name=args.name
    )

    user_args = {}
    if args.o:
        for kwarg in args.o:
            kw, arg = kwarg.split(":")
            user_args[kw] = arg

    api.set_qsub_args(jobname=args.name, **user_args)
    jobid = api.add_to_que()
    print(f"Added to que: {api.name} / {jobid}")





    # jobname=jobname,
    # output = os.path.join("logs", jobname),
    # error = os.path.join("logs",jobname)