# Screener
Repository for Master thesis "K-mer based quantification of biosynthetic gene clusters"

This repository holds the used in the project, as project work is more fluid so is the structure of the repository.
The code is focused on processing and analysing the output generated by the CAMISIM analysis tool and is therefore,
less optimized for analyzing outside samples. 
This is a known flaw with the code and will be rectified.

## Disclaimer about environment.

The code was developed and run on the Computerome health-tech supercomputer.
The implementations are therefore dependent on the que-system and module structure of the supercomputer.
This is especially true for runnning the full pipelines.

Scripts and the python classes (lib/*) can be used independent of the Computerome infrastructure.

## Running a pipeline.

Assuming working on Computerome.

`
git clone git@github.com:HenrikSpiegel/Screener.git
cd Screener 
module purge
module load tools anaconda3/2021.11       # a newer python is required to run the pipeline.
python -m pipeline.data_simulation_small  # -d flag (dryrun) will generate dependency/status graph in log/pipeline without adding jobs
`

Note for running times.
The development set (pipeline.data_simulation_small) finishes within hours
The benchmark set (pipeline.data_simulation_large) takes 2-3 days to run dependending on traffic.

The CAMISIM simulation steps are generally the slowest.


## Repository structure

The project is structured into multiple submodules.

In order to comply with the computerome quesystem, the majority code was written as scripts with command line interface (CLI)
The qsub-modules and pipelines are an exception.

### analysis:
Specific scripts for creating various analysis (primarily visual). 

### config
Config files for the full pipelines.

### data:
sugggested data folder - ships empty.

### lib:
High level python modules with large classes handling the initial catalogue assembly (catalogue_assembler.py) and later refinement (magpy.py)

Note: It was found that the catalouge_assembler scaled poorly with large number of BGC clusters.
A new method which is mostly based on the Jellyfish2 infrastructure was implemented instead which is much faster,
and scales better as parallelization is easily achieved.

The Jellyfish implementation is included in scripts/generate_catalogue_jellyfish.py, however, it was not used
in the main pipelines as the problem was only observed for a larger external dataset, and the python implementation
was suitable for the smaller sets.

### logs:
log output from pipelines and individual qsub (computerome que-system) jobs.

### Pipeline:
pipeline_base.py includes a base pipeline class which can processes jobs with dependencies and provide graphical pipeline progressions.
The pipeline_base class is heavily linked with the qsub_modules and rely on the qsub_module classes to report job-status etc.

The individual pipelines can the be created by adding jobs (qsub_modules) to the pipeline and declaring their dependencies.

data_simulation_large.py is the full pipeline described in the master thesis as the benchmark set, and the data_simulation_small.py describes the creation
and analysis of the development set.

While the pipelines can be run as a qsub-job themselves (via. qsub_modules.add_to_que) it is recommoneded to run them interactively (preferable in a screen) as logging
output might be relevant. 

The pipeline can be restarted at a breakpoint as individual jobs are marked complete with `.success` files.
Likewise, if a step is needed to be rerun -> the relevant .success file can be deleted.

The maximum number of simultaniously running jobs can be adjusted when setting up the pipeline.


### qsub_modules
qsub_modules.base contains the base class which have methods for adding qsub job-scripts to the computerome que and monitor their status.  

The qsub_modules generally wrap a script (scripts/..) and includeds a preflight method for input checking aswell as the qsub methods.
For some qsub_modules the job is defined in the module and not as a stand-alone script.  
This is generally the case when calling other tools, such as simple calls to samtools 

qsub_modules.antismash is an example a job, the class inherits from the qsub base class and defines required modules, resources and the actual commands.
The module can be run as a standalone, but is developed to be implemented in a pipeline.

### scripts

Various low-level scripts both including stand-alone scripts as well as calling various other tools.


```
.
├── analysis
│   ├── 01_compare_input_bgc_genera.py
│   ├── 01_compare_input_bgc.py
│   ├── 02_plot_bgc_coverage_simulations.py
│   ├── ...
├── analysis_notebooks
│   └── 11_mag_test_subset.ipynb
├── config
│   ├── analysis_ibdmdb.ini
│   ├── project_config_init.ini
│   └── project_config_large.ini
├── data
│   ├── simulated_data_init
│   ├── ...
├── lib
│   ├── catalogue_assembler.py
│   ├── magpy.py
│   └── utils.py
├── logs
│   ├── pipeline
│   └── qsub
├── pipeline
│   ├── data_analysis.py
│   ├── data_simulation_large.py
│   ├── data_simulation_small.py
│   ├── pipeline_base.py
├── qsub_modules
│   ├── add_to_que.py
│   ├── antismash.py
│   ├── base.py
│   ├── blastn_pw.py
│   ├── camisim.py
│   ├── ...
├── scripts
│   ├── antismash_as_fasta.py
│   ├── average_readlength.py
│   ├── blast_demo_prep.py
│   ├── camisim_combine_descriptions.py
│   ├── camisim_describe_run.py
│   ├── collect_count_matrices.py
│   ├── ...

```



