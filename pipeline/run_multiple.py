import os, sys
import subprocess
import time
import numpy as np
gbs_to_run = np.round([x*10**-y for x in range(1,10) for y in (1,2,3)] + [1], 3) #Note rounding should equal max powering
gbs_to_run.sort()

id_map = dict()
for gb in gbs_to_run:
    print("--"+str(gb), file=sys.stderr)
    batcmd = f"python -m pipeline.simulate_x_gb --readsGB {gb}"
    result = subprocess.check_output(batcmd, shell=True)
    id_map[str(gb)]=result
    time.sleep(1)

print(id_map)