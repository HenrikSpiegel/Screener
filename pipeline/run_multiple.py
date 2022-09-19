import os
import time
gbs_to_run = [x*10**-y for x in (1,2) for y in (0,1,2,3)]

for gb in gbs_to_run:
    os.system(f"python -m pipeline.simulate_x_gb --readsGB {gb}")
    time.sleep(1)