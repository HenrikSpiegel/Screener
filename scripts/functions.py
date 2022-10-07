import subprocess
import sys, os
import time
from typing import Callable

def syscall(command):
    job = subprocess.run(command.split(),
          stdout=subprocess.PIPE, universal_newlines=True)
    result = job.stdout.split('\n')
    return result

def submit2(command, runtime, cores, ram, directory='', modules='', group='dtu_00009',
    output='/dev/null', error='/dev/null', jobname="henspi", email='', dependency = [], test=False):
    """
    Function to submit a job to the Queueing System - without jobscript file
    Parameters are:
    command:   The command/program you want executed together with any parameters.
               Must use full path unless the directory is given and program is there.
    directory: Working directory - where should your program run, place of your data.
               If not specified, uses current directory.
    modules:   String of space separated modules needed for the run.
    runtime:   Time in minutes set aside for execution of the job.
    cores:     How many cores are used for the job.
    ram:       How much memory in GB is used for the job.
    group:     Accounting - which group pays for the compute.
    output:    Output file of your job.
    error:     Error file of your job.
    dependency: list of qsub ids which must be ok before it can run.
    """
    runtime = int(runtime)
    cores = int(cores)
    ram = int(ram)
    if not isinstance(dependency, list):
        dependency = [dependency]
    if cores > 38:
        print("Can't use more than 38 cores on a node")
        sys.exit(1)
    if ram > 188:
        print("Can't use more than 120 GB on a node")
        sys.exit(1)
    if runtime < 1:
        print("Must allocate at least 1 minute runtime")
        sys.exit(1)
    minutes = runtime % 60
    hours = int(runtime/60)
    walltime = "{:d}:{:02d}:00".format(hours, minutes)
    if directory == '':
        directory = os.getcwd()
    # Making a jobscript
    script = '#!/bin/sh\n'
    script += '#PBS -A ' + group + ' -W group_list=' + group + '\n'
    script += '#PBS -e ' + error + ' -o ' + output + '\n'
    script += '#PBS -d ' + directory + '\n'
    script += '#PBS -l nodes=1:ppn=' + str(cores) + ',mem=' + str(ram) + 'GB' + '\n'
    script += '#PBS -l walltime=' + walltime + '\n'
    script += '#PBS -N '+ jobname + '\n'
    if dependency != []:
        script += '#PBS -W depend=afterok:'+':'.join([str(id) for id in dependency]) + '\n'
    if email != '':
        script += '#PBS -M '+ email + '\n'
        script += '#PBS -m ae' + '\n'
    script += "set -e" +'\n' #if any calls fail the whole script return non-zero.
    if modules != '':
        script += 'module purge' + '\n'
        script += 'module load ' + modules + '\n'
    script += "set -e" +'\n' #if any calls fail the whole script return non-zero.
    script += command + '\n'
    # The submit
    if test:
        print(script)
        return("NoID")
    else:
        job = subprocess.run(['qsub'], input=script, stdout=subprocess.PIPE, universal_newlines=True)
        jobid = job.stdout.split('.')[0]
        return jobid.strip()

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


def check_finish(jobid: int, _success_func: Callable, wait_for_finish=False, wait_time=120, **kwargs_success_func):
    # If the job is finished sucessfully the qsub status will be C
    # AND the stderr will contain "Metagenome simulation finished"
    # By convention the submit2() the error output will be located in {wdir}/logs/jobname_stderr.txt
    # The jobname can be found from the qsub output
    import subprocess
    

    #loop waiting for finished job:
    job_is_finished = job_finished(jobid)
    while not job_is_finished:
        print("Job not finished", file=sys.stderr)

        if not wait_for_finish:
            print("Not waiting for finish", file=sys.stderr)
            return 2
        print(f"sleeping ({wait_time}s)", file=sys.stderr)
        time.sleep(wait_time)
        job_is_finished = job_finished(jobid)

    if _success_func(**kwargs_success_func):
        print("sucessfull run", file=sys.stderr)
        return 0
    if not _success_func(**kwargs_success_func):
        print("Failed run", file=sys.stderr)
        return 1