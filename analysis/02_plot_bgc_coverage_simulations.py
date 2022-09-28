#!/usr/bin/env python3
import os
import glob
from typing import List
from Bio import SeqIO
import numpy as np
import pathlib
import pandas as pd
from tqdm import tqdm
import plotly.express as px


def gen_rrna_annot(barrnap_file: str = "results/man_rrna_detection_all_bgc.txt") -> dict:
    annot = dict()
    for line in open(barrnap_file, "r"):
        if line.startswith("#"):
            continue
        line_seq = line.strip().split("\t")
        feature_dict = {
            'type': 'rRNA',
            'location': (int(line_seq[3]), int(line_seq[4])),
            'locus_tag': line_seq[8].split("=")[1].split(";")[0]
        }
        entry = annot.get(line_seq[0], [])
        entry.append(feature_dict)
        annot[line_seq[0]] = entry
    return annot 

def gen_pileup_fig(df, dataset_label="testData", bgc="test_bgc", expected_cov:float = None):
    fig = px.area(df, x="pos",y="cov",
                  title=f"Coverage plot for {bgc} <br><sup>{dataset_label} dataset</sup>",
                  labels={
                     "pos": "Position in Cluster",
                     "cov": "Coverage (count)"
                 },
                 template="plotly_white"
                 )
    if expected_cov:
        fig.add_hline(y=expected_cov, fillcolor="red")
    return fig

def add_genomic_annotation(fig, genbank_record, custom_annot: List[dict]=[]):
    y1_relative = 0.12
    anno_yshift_relative = 0.85
    color_map = {
        "CDS":"grey",
        "rRNA": "darkred"
    }

    y_magnitude = fig.data[0]["y"].max()
    shape_y1 = -(y1_relative * y_magnitude)
    anno_y1  = -(y1_relative * y_magnitude)
    anno_yshift = -(anno_yshift_relative*y_magnitude)

    for feat in genbank_record.features:
        if feat.type=="region":
            pass
        if feat.type=="CDS":
            location = (min(feat.location), max(feat.location))
            locus_tag = feat.qualifiers["locus_tag"]
            fig.add_shape(
                name=locus_tag[0],
                x0 = location[0],
                x1 = location[1],
                y0=-0.1,
                y1=shape_y1,
                yref="y",
                fillcolor=color_map[feat.type]
            )
            fig.add_annotation(
                text = locus_tag[0],
                x = np.mean(location),
                y=anno_y1,
                yshift=-35,
                showarrow=False,
                textangle=90,
                align="right"
            )
    for entry in custom_annot: #does nothing if custom_annot == []
        fig.add_shape(
                name=entry["locus_tag"],
                x0 = entry["location"][0],
                x1 = entry["location"][1],
                y0=-0.1,
                y1=shape_y1,
                yref="y",
                fillcolor=color_map[entry["type"]]
        )
        fig.add_annotation(
            text = entry["locus_tag"],
            x = np.mean(entry["location"]),
            y=anno_y1,
            yshift=-35,
            showarrow=False,
            textangle=90,
            align="right"
        )
    return fig


if __name__ == "__main__":
    import sys
    sys.stderr.write("--- Running: "+__file__+"\n")

    outdir = "results/02_bgc_coverage"
    coverage_files = "data/simulated_data/quantification_map/*GB/coverage.tsv"
    genbank_dir = "data/simulated_data/antismash/input_genomes"

    custom_annotation = gen_rrna_annot("results/man_rrna_detection_all_bgc.txt")

    # For each dataset (readsGB)
    for fp in tqdm(glob.glob(coverage_files)):
        label = os.path.basename(os.path.dirname(fp))
        pathlib.Path(os.path.join(outdir, label)).mkdir(parents=True, exist_ok=True)

        df =  pd.read_csv(fp, sep="\t", names=["bgc","pos","cov"])
        df_grouped = df.groupby("bgc")
        # For each genome:bgc entry in dataset
        for name, group in tqdm(df_grouped, leave=False):
            outfile = os.path.join(outdir, label, name+".png")
            df = group.copy().sort_values("pos").reset_index(drop=True)
            fig = gen_pileup_fig(df, dataset_label=label, bgc=name.replace("_","."))
            
            #get annotation file:
            genbank_fp = os.path.join(genbank_dir, f"{name}.gbk")
            record = next(SeqIO.parse(genbank_fp, "genbank")) #Assumes only one record

            #check for custome annotations
            custom_entry = custom_annotation.get(name, [])
            
            #Add annotations a write file.
            fig_anno = add_genomic_annotation(fig, record, custom_entry)
            fig.write_image(outfile, height=700, width=1000)
