import os
import subprocess
import time
gbs_to_run = [x*10**-y for x in range(1,10) for y in (1,2,3)] + [1]

id_map = dict()
for gb in gbs_to_run:
    batcmd = f"python -m pipeline.simulate_x_gb --readsGB {gb}"
    result = subprocess.check_output(batcmd, shell=True)
    id_map[str(gb)]=result
    time.sleep(1)

print(id_map)