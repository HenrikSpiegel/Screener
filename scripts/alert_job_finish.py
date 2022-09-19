import subprocess
def job_finished(jobid):
    def _get_qstat(jobid) -> dict:
        call_qstat = f"qstat -f {jobid}"
        output_stream = subprocess.run(call_qstat.split(" "), stdout=subprocess.PIPE).stdout.decode('utf-8')
        output_cleaned = [x.strip() for x in output_stream.replace("\n\t", "").split("\n")]

        output_dict = {}
        for line in output_cleaned[1::]:
            if not line or "=" not in line:
                continue
            k, v = line.split(" = ")
            if k in ["Error_Path", "Output_Path"]:
                v = v.split(":")[1]
            output_dict[k] = v
        return output_dict
    qstat = _get_qstat(jobid)
    return qstat['job_state'] == "C"

if __name__ == "__main__":
    import sys
    import time
    id = sys.argv[1]

    while not job_finished(id):
        time.sleep(10)
    print(f"Job finished: "+str(id), file=sys.stderr)