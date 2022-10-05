import os, sys
import subprocess
import time
import numpy as np
import configparser

config = configparser.ConfigParser()
config.read("config/project_config.ini")

if config.get("Simulation", "ReadsGBStep") != "":
    rstep = config.getfloat("Simulation", "ReadsGBStep")
    rmin = config.getfloat("Simulation", "ReadsGBMin")
    rmax = config.getfloat("Simulation", "ReadsGBMax")
    precision = config.get("Simulation", "ReadsGBStep").split(".")[1].__len__()
    gbs_to_run = np.round(np.arange(rmin, rmax+rstep, rstep), precision).tolist()
else:
    gbs_to_run = []

if config.get("Simulation", "ReadsGBExtra", fallback=None):
    gbs_extra = [float(x.strip()) for x in config.get("Simulation", "ReadsGBExtra").split()]
    gbs_to_run.extend(gbs_extra)


if __name__ == "__main__":
    id_map = dict()
    for gb in gbs_to_run:
        print("--"+str(gb), file=sys.stderr)
        batcmd = f"python -m pipeline.simulate_x_gb --readsGB {gb}"
        result = subprocess.check_output(batcmd)
        id_map[str(gb)]=result.decode().strip()