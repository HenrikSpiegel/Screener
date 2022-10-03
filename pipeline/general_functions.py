import logging
import pathlib
import sys, os
import subprocess

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

    if not isinstance(dependency, list):
        dependency = [dependency]
    runtime = int(runtime)
    cores = int(cores)
    ram = int(ram)
    if cores > 38:
        print("Can't use more than 30 cores on a node")
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
    if modules != '':
        script += 'module load ' + modules + '\n'

    script += command + '\n'
    # The submit
    if test:
        print(script)
        return("NoID")
    else:
        job = subprocess.run(['qsub'], input=script, stdout=subprocess.PIPE, universal_newlines=True)
        jobid = job.stdout.split('.')[0]
        return jobid.strip()


def get_log(jobtag="test", log_name=None, logdir = "logs/pipeline", lvl=logging.DEBUG):
    if not log_name:
        log_name = os.path.basename(__file__).split(".")[0]+"_"+jobtag
    logger = logging.getLogger(log_name)
    
    #remove old handlers:
    while logger.handlers:
        logger.removeHandler(logger.handlers[0])
    
    logger.setLevel(lvl)
    F = "[%(asctime)s %(name)s:%(funcName)s]%(levelname)s: %(message)s"
    formatter = logging.Formatter(F, datefmt='%d-%b-%y %H:%M:%S')
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    #stream_handler.setLevel(lvl)
    logger.addHandler(stream_handler)
    
    if logdir:
        fp_log = os.path.join(logdir, log_name)
        if os.path.isfile(fp_log): 
            os.unlink(fp_log)
            
        pathlib.Path(os.path.dirname(fp_log)).mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(fp_log)
        file_handler.setFormatter(formatter)
        #file_handler.setLevel(lvl)
        logger.addHandler(file_handler)
        logger.debug("logfile at -> "+fp_log)
    return logger

def raise_w_log(log, exception:Exception, msg):
    log.error(msg)
    raise exception(msg)