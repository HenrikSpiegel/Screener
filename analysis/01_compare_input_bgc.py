
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
    fp_pw_blast  = Path("data/simulated_data/blast_pairwise/input_bgc/combined_blast_results.tsv")
    output       = Path("results/") / Path(__file__).stem
    output.mkdir(parents=True, exist_ok=True)

    assert(antismashdir.exists())
    assert(fp_pw_blast.exists())

    gbk_files = [x for x in antismashdir.glob("*.gbk") if x.name != "combined.gbk"]
    assert(len(gbk_files) > 0)

    seq_annotation = annotation_from_genbank(gbk_files)

    df_blast = pd.read_csv(fp_pw_blast, sep="\t")
  

    #Get summarised p_ident of query coverage.
    df_blast["weighted_pident"] = df_blast.pident * df_blast.length
    df_grouped = df_blast.groupby(["qaccver", "saccver"])
    summarised = []
    for name, group in df_blast.groupby(["qaccver", "saccver"]):
        qcov = group.qcovs.iloc[0]
        summarised_pident = group.weighted_pident.sum() / group.length.sum()
        summarised.append(name + (qcov, summarised_pident))
    df_pairwise = pd.DataFrame(
        data=summarised,
        columns = ("query", "subject", "query_coverage", "summarised_pident")
    )

    df_meta = df_pairwise.assign(
        query_type = [seq_annotation[x]["assigned_type"] for x in df_pairwise['query']],
        subject_type = [seq_annotation[x]["assigned_type"] for x in df_pairwise['subject']],
        query_len  = [seq_annotation[x]["seq_len"] for x in df_pairwise['query']],
        subject_len  = [seq_annotation[x]["seq_len"] for x in df_pairwise['subject']]
    )
    
    #Calculate pairwise similarity as well as aligned bases.
    ### We have the following similarity measure:

    #sim(q, s) = qcov*qlen *sum_pident / max(qlen, sublen)

    #In essense we find the ratio between aligned bases (scaled by the % identity of aligned bases.)
    #and the maximum possible aligned bases
    def similarity_metric(row):
        aligned_bases = row['query_coverage']*row['query_len']*(row['summarised_pident']/100)
        maximum_aligned_bases =   row[['query_len', 'subject_len']].max()
        return (aligned_bases/(1*10**5), aligned_bases / maximum_aligned_bases)

    df_meta[["aligned_bases", "similarity"]] = df_meta.apply(similarity_metric, axis=1, result_type="expand")


    ## Visualization
    # get value matrix:
    df_wide = df_meta.pivot_table(index="query", columns="subject", values="similarity")

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
        height=850
    )
    fig.update_xaxes(tickangle=90, ticklabeloverflow="allow")

    print(f"Writing outputs to : {output}", file=sys.stderr)
    fig.write_image(file=output/'blast_bgc_identity_heatmap.png', height=850, width=1000)
    df_meta.to_csv(output/'blast_annotated_pw_bgc_comparison.csv')

    print("Finished: " + __file__, file=sys.stderr)
    

