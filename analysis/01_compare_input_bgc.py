#!/usr/bin/env python3

# This script runs comparative analysis on a set of bgcs.
# Output is a csv with pairwise identities.
import os
import re
import sys
import time
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from typing import List
from scripts.functions import submit2, check_finish,job_finished
import scripts.qsub_needleman as needleman
from Bio import SeqIO

#workingdir = "/home/projects/dtu_00009/people/henspi/git/AntibioticaScreening/project"

#### Find relevant files.
antismashdir = "data/simulated_data/antismash/input_genomes"
gbk_files = [x for x in os.listdir(antismashdir) if x.endswith(".gbk") and x != "combined.gbk"]

#### run pairwise alignment:
# we want to run all the pairwise combinations. 

combinations = [(gbk_files[i], gbk_files[i+1::]) for i in range(len(gbk_files)-1)] #if gbk_files=[1,2,3] ->  [(1, [2,3]), (2,[3])]

combined_call = []
output_fhs = []
for comb in combinations:
    seq1 = comb[0]
    for seq2 in comb[1]:
        output_fn = seq1.replace(".gbk","")+"_"+seq2.replace(".gbk",".align")
        output_fh = os.path.join("data/simulated_data/pw_needle", output_fn)
        output_fhs.append(output_fh)
        seq1_fh = os.path.join(antismashdir, seq1)
        seq2_fh = os.path.join(antismashdir, seq2)
        needleman.preflight(seq1_fp=seq1_fh, seq2_fp=seq2_fh, outfile=output_fh)
        call = needleman.generate_syscall(seq1_fp=seq1_fh, seq2_fp=seq2_fh, outfile=output_fh)
        combined_call.append(call)

#check if we already have the output files
output_exists = [os.path.exists(p) for p in output_fhs]
if any(output_exists):
    print(f"Found already run outputs ({sum(output_exists)}), remaining outputs ({len(output_exists)-sum(output_exists)})", file=sys.stderr)

    combined_call = [call for call, exist in zip(combined_call, output_exists) if not exist]

#split calls.
calls_per_job = 10
combined_call_split = [combined_call[x:x+calls_per_job] for x in range(0, len(combined_call), calls_per_job)]

print(f"Running comp_bgc as ({len(combined_call_split)}) subjobs")
runids = []
for i, subcall in enumerate(combined_call_split):
    runid = submit2(
        command="\n".join(subcall),
        **needleman.gen_qsub_args(job_tag=str(i))
    )
    print(f"compare_bgs runid: {runid}", file=sys.stderr)
    runids.append(runid)
    time.sleep(1) #Lets be nice to the que system.

if runids != []:
    # Wait for jobs to finish.
    job_has_finished = [job_finished(x) for x in runids]
    while not all(job_has_finished):
        jobs_not_finished = [job for job, is_finished in zip(runids, job_has_finished) if not is_finished]
        print(f"waiting for ({len(jobs_not_finished)}) jobs", file=sys.stderr)
        time.sleep(20)
        job_has_finished = [job_finished(x) for x in runids]

#check all output files exists
def is_success(output_expected: List[str]) -> bool:
    output_files_found = [os.path.isfile(x) for x in output_expected]
    print(f"is_success: {sum(output_files_found)}/{len(output_files_found)} files found", file=sys.stderr)
    return (all([os.path.isfile(x) for x in output_expected]), list(set(output_expected).difference(set(output_files_found))))

job_succesful, missing_files = is_success(output_fhs)

if not job_succesful:
    raise RuntimeError(f"Failed Getting all outputs: missing files -> {missing_files}")

############################### DATA GRABBING AND PRESENTING ###############################

## Gather results and publish.

### Grab metadata from antismash genbank files:
seq_annotation = {}
for file in gbk_files:
    seqname = file[:-4]
    file_fp = os.path.join(antismashdir, file)
    
    record = next(SeqIO.parse(file_fp, "genbank")) #assume 1 record per file.
    seq = record.seq.__str__()
    seq_len = len(record.seq)
    #grab predicted product.
    for feat in record.features:
        if feat.type == "protocluster":
            assigned = feat.qualifiers["product"]
            break
    seq_annotation[seqname] = {
        "fp": file_fp,
        "seq": seq,
        "seq_len": seq_len,
        "assigned_type": assigned
    }

### Grab pairwise alignment files.
results = []
for fn in output_fhs:
    for line in open(fn, "r"):
        if line.startswith("#    -asequence"):
            seq1 = line.strip().rsplit("/",1)[1][:-4]
            continue
        elif line.startswith("#    -bsequence"):
            seq2 = line.strip().rsplit("/",1)[1][:-4]
            continue
        elif line.startswith("# Identity:"):
            reg_res = re.search(r"\((.+)%\)", line)
            if reg_res:
                iden = float(reg_res.group(1))
            break
    results.append((seq1, seq2, iden))

df = pd.DataFrame(columns = ["seq1","seq2","iden"], data=results)

### Combine sources to create master df

# duplicates to go from triangle to full pw matrix
df_combined = pd.concat([
    pd.DataFrame(columns = ["seq1","seq2","iden"], data=results),
    pd.DataFrame(columns = ["seq2","seq1","iden"], data=results) 
]) 

## Add metadata
df_meta = df_combined.assign(
    seq1_type = [seq_annotation[x]["assigned_type"] for x in df_combined.seq1],
    seq2_type = [seq_annotation[x]["assigned_type"] for x in df_combined.seq2],
    seq1_len  = [seq_annotation[x]["seq_len"] for x in df_combined.seq1],
    seq2_len  = [seq_annotation[x]["seq_len"] for x in df_combined.seq2]
)
df_fp = "results/01_annotated_pw_bgc_comparison.csv"
print("Data written to -> "+df_fp)
df_meta.to_csv(df_fp)

### plotting
df_wide = df_combined.pivot_table(index="seq1", columns="seq2", values="iden")
x_ticks = df_wide.columns
y_ticks = df_wide.index
values = df_wide.values
# Generate layered x and y ticks # TODO: fix outer x-ticks overlapping.

refseq2species = {
    "NZ_LT906445.1": "V. Parvula",
    "NZ_CP053893.1": "C. beijerinckii",
    "NC_014328.1":   "C. ljungdahlii",
    "NZ_CP020566.1": "V. atypica",
    "NZ_LT906470.1":   "V. rodentium "
}

x_ticks_expanded = []
y_ticks_expanded = []
test = 0
for x, y in zip(x_ticks, y_ticks):
    genome_x, region_x = x.rsplit(".",1)
    genome_x = refseq2species[genome_x]
    genome_y, region_y = y.rsplit(".",1)
    genome_y = refseq2species[genome_y]
    type_x = "; ".join(seq_annotation[x]["assigned_type"]) + "_Reg"+ region_x[-3:] + f"_{seq_annotation[y]['seq_len']/1000:.1f}kb"
    type_y = "; ".join(seq_annotation[y]["assigned_type"]) + "_Reg"+ region_y[-3:] + f"_{seq_annotation[y]['seq_len']/1000:.1f}kb"
    x_ticks_expanded.append([genome_x, type_x])
    y_ticks_expanded.append([genome_y, type_y])
    test += 1

x_ticks_expanded = np.array(x_ticks_expanded).T
y_ticks_expanded = np.array(y_ticks_expanded).T 

fig = go.Figure(data=go.Heatmap(
    z=values,
    x=x_ticks_expanded,
    y=y_ticks_expanded,
    colorscale="Viridis",
    colorbar=dict(title='Seq identity (%)')
))

fig.update_layout(
    height=850
)
fig.update_xaxes(tickangle=90, ticklabeloverflow="allow")
fig_fh = "results/01_bgc_heatmap.png"
print("Heatmap written to -> "+fig_fh)
fig.write_image(fig_fh, height=850, width=1000)


# TODO: check_finish needs to accept multiple ids or perhaps have a better syntax for child processes.
#finished_successful = check_finish(runid, _success_func=is_success, wait_for_finish=True, wait_time=30, output_files=output_fhs)
#print(finished_successful, file = sys.stderr)