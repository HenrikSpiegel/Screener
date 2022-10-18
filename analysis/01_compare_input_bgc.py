
import sys
from Bio import SeqIO
import pandas as pd
from typing import List
from pathlib import Path

import numpy as np
import plotly.graph_objects as go

def annotation_from_genbank(gbk_files: List[Path]) -> dict:
    ### get seqlength for input files.
    seq_annotation = {}
    for file in gbk_files:
        seqname = file.stem

        record = next(SeqIO.parse(file, "genbank")) #assume 1 record per file.
        seq = record.seq.__str__()
        seq_len = len(record.seq)
        #grab predicted product.
        for feat in record.features:
            if feat.type == "protocluster":
                assigned = feat.qualifiers["product"]
                break
        seq_annotation[seqname] = {
            "fp": file,
            "seq": seq,
            "seq_len": seq_len,
            "assigned_type": assigned
        }
    return seq_annotation

if __name__ == '__main__':

    print("Running: " + __file__, file=sys.stderr)

    antismashdir = Path("data/simulated_data/antismash/input_genomes")
    fp_pw_blast  = Path("data/simulated_data/pairwise/combined_blast_results.tsv")
    output       = Path("results/") / Path(__file__).stem
    output.mkdir(parents=True, exist_ok=True)

    assert(antismashdir.exists())
    assert(fp_pw_blast.exists())

    gbk_files = [x for x in antismashdir.glob("*.gbk") if x.name != "combined.gbk"]
    assert(len(gbk_files) > 0)

    seq_annotation = annotation_from_genbank(gbk_files)

    df_pairwise = pd.read_csv(fp_pw_blast, sep="\t")
    df_pairwise.rename(columns={'qaccver':'seq1', 'saccver':'seq2'}, inplace=True)

    df_meta = df_pairwise.assign(
    seq1_type = [seq_annotation[x]["assigned_type"] for x in df_pairwise.seq1],
    seq2_type = [seq_annotation[x]["assigned_type"] for x in df_pairwise.seq2],
    seq1_len  = [seq_annotation[x]["seq_len"] for x in df_pairwise.seq1],
    seq2_len  = [seq_annotation[x]["seq_len"] for x in df_pairwise.seq2]
    )
    df_meta

    # get value matrix:
    df_wide = df_pairwise.pivot_table(index="seq1", columns="seq2", values="pident")
    x_ticks = df_wide.columns
    y_ticks = df_wide.index
    values = df_wide.values

    # Generate nested ticks.
    x_ticks_expanded = []
    y_ticks_expanded = []
    test = 0
    for x, y in zip(x_ticks, y_ticks):
        genome_x, region_x = x.rsplit(".",1)
        genome_y, region_y = y.rsplit(".",1)
        type_x = "; ".join(seq_annotation[x]["assigned_type"]) + "_Reg"+ region_x[-3:] + f"_{seq_annotation[y]['seq_len']/1000:.1f}kb"
        type_y = "; ".join(seq_annotation[y]["assigned_type"]) + "_Reg"+ region_y[-3:] + f"_{seq_annotation[y]['seq_len']/1000:.1f}kb"
        x_ticks_expanded.append([genome_x, type_x])
        y_ticks_expanded.append([genome_y, type_y])
        test += 1

    x_ticks_expanded = np.array(x_ticks_expanded).T
    y_ticks_expanded = np.array(y_ticks_expanded).T 

    # Generate plot
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

    print(f"Writing outputs to : {output}", file=sys.stderr)
    fig.write_image(file=output/'bgc_identity_heatmap.png', height=850, width=1000)
    df_meta.to_csv(output/'annotated_pw_bgc_comparison.csv')

    print("Finished: " + __file__, file=sys.stderr)
    

