import os, sys
import subprocess
import time
import numpy as np
import configparser
import json

config = configparser.ConfigParser()
config.read("config/project_config.ini")

if config.get("Simulation", "ReadsGBStep") != "":
    rstep = config.getfloat("Simulation", "ReadsGBStep")
    rmin = config.getfloat("Simulation", "ReadsGBMin")
    rmax = config.getfloat("Simulation", "ReadsGBMax")
    precision = config.get("Simulation", "ReadsGBStep").split(".")[1].__len__() #Høker significant digits.
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
        sim_cmd = f"python -m pipeline.run_simulation --readsGB {gb}"
        sim_out = subprocess.check_output(sim_cmd.split(" "))
        gb_ids = json.loads(sim_out.decode().strip())

        eval_cmd_list = ["python", "-m", "pipeline.run_quantification_sim",
        "--readsGB", str(gb),
        "--dependencies", json.dumps(gb_ids),
        "--runQuantifierKmer", 
        "--runQuantifierMap"]
        #eval_cmd = f"python -m pipeline.run_quantification_sim --readsGB {gb} --dependencies '{json.dumps(gb_ids)}' --runQuantfierKmer --runQuantifierMap "
        try:
            eval_out = subprocess.check_output(eval_cmd_list)#, shell=True)
        except Exception as err:
            print(err)
            print(" ".join(eval_cmd_list))
            raise err


        eval_ids = json.loads(eval_out.decode().strip())
        gb_ids.update(eval_ids)

        id_map[str(gb)]=gb_ids
    print(json.dumps(id_map))