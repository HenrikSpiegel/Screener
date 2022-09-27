import os, sys
import subprocess
import time
import numpy as np
import configparser

config = configparser.ConfigParser()
config.read("config/project_config.ini")

rstep = config.getfloat("Simulation", "ReadsGBStep")
rmin = config.getfloat("Simulation", "ReadsGBMin")
rmax = config.getfloat("Simulation", "ReadsGBMax")
precision = config.get("Simulation", "ReadsGBStep").split(".")[1].__len__()
gbs_to_run = np.round(np.arange(rmin, rmax+rstep, rstep), precision).tolist()

if config.get("Simulation", "ReadsGBExtra", fallback=None):
    gbs_extra = [float(x) for x in config.get("Simulation", "ReadsGBExtra").split(",")]
    gbs_to_run.extend(gbs_extra)


if __name__ == "__main__":
    id_map = dict()
    for gb in gbs_to_run:
        print("--"+str(gb), file=sys.stderr)
        batcmd = f"python -m pipeline.simulate_x_gb --readsGB {gb} --runQuantifierMap --runQuantifierKmer"
        result = subprocess.check_output(batcmd, shell=True)
        id_map[str(gb)]=result.decode().strip()
        time.sleep(1)

    sys.stdout.write(str(id_map)) #TODO: Use json dump instead.