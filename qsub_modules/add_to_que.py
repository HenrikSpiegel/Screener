import argparse
from qsub_modules.base import Base
from pathlib import Path
import configparser
import logging




class AddToQue(Base):
    def __init__(self, command:str, success_file:Path, name:str=None, loglvl = "DEBUG"):
        self.command = command
        self.name = name or "AddToQue"
        self.jobtag = name
        self.loglvl = loglvl

        if success_file:
            self.success_file = Path(success_file)

    qsub_requirements = dict(
        modules = "tools anaconda3/2021.11",
        runtime = 360,
        cores = 10,
        ram=50,
        group="dtu_00009",
    )

    @property
    def log_setup(self):
        return dict(
                name = self.name,
                level = self.loglvl,
                log_file = None
            )    

    def generate_syscall(self):
        call = self.command
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
        #test = args.test,
        name=args.name
    )

    user_args = {}
    if args.options:
        for kwarg in args.options:
            kw, arg = kwarg.split(":")
            user_args[kw] = arg

    api.set_qsub_args(jobname=args.name, **user_args)
    api.add_to_que(test=args.test)
