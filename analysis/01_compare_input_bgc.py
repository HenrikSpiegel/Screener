
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

    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument("--antismashdir", required=True, type =Path)
    parser.add_argument("--blast-table", required=True, type=Path)
    parser.add_argument("-o", required=True, type=Path)
    args = parser.parse_args()

    antismashdir = args.antismashdir #Path("data/simulated_data/antismash/input_genomes")
    symmetric_blastres = args.blast_table #Path("data/simulated_data/blast_pairwise/input_bgc/pairwise_table_symmetric.tsv")
    output = Path(args.o) #Path("results") / Path(__file__).stem
    output.mkdir(parents=True, exist_ok=True)
    
    gbk_files = [x for x in antismashdir.glob("*.gbk") if x.name != "combined.gbk"]
    seq_annotation = annotation_from_genbank(gbk_files)

    df_wide = pd.read_csv(symmetric_blastres, sep="\t", index_col=0)
    
    x_ticks = df_wide.columns
    y_ticks = df_wide.index
    values_similarity = df_wide.round(2).fillna("").values

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
        z=values_similarity,
        text=values_similarity,
        texttemplate="%{text}",
        textfont={"size":10},

        x=x_ticks_expanded,
        y=y_ticks_expanded,
        colorscale="Inferno",
        colorbar=dict(title='Similarity (%)'
                    )
    ))

    fig.update_layout(
        height=850,
        title = "Similarity between BGCs of the five input genomes<br><sup>Similarity ratio of perfectly aligned bases to maximum possible aligned bases."
    )
    fig.update_xaxes(tickangle=90, ticklabeloverflow="allow")

    print(f"Writing outputs to : {output}", file=sys.stderr)
    fig.write_image(file=output/'blast_bgc_identity_heatmap.png', height=850, width=1000)

    print("Finished: " + __file__, file=sys.stderr)
    

