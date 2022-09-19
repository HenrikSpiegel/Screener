#!/usr/bin/env python3
import os
import glob
from Bio import SeqIO
import pathlib
import pandas as pd
from tqdm import tqdm
import plotly.express as px

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

def add_genomic_annotation(fig, genbank_record):
    y1_relative = 0.12
    anno_yshift_relative = 0.85

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
                fillcolor="grey"
            )
            fig.add_annotation(
                text = locus_tag[0],
                x = sum(location)/2,
                y=anno_y1,
                yshift=-35,
                showarrow=False,
                textangle=90,
                align="right"
            )
    return fig


if __name__ == "__main__":

    outdir = "results/02_bgc_coverage"
    coverage_files = "data/simulated_data/quantification/*GB/coverage.tsv"
    genbank_dir = "data/simulated_data/antismash/input_genomes"

    for fp in tqdm(glob.glob(coverage_files)):
        label = os.path.basename(os.path.dirname(fp))
        pathlib.Path(os.path.join(outdir, label)).mkdir(parents=True, exist_ok=True)

        df =  pd.read_csv(fp, sep="\t", names=["bgc","pos","cov"])
        df_grouped = df.groupby("bgc")
        for name, group in tqdm(df_grouped, leave=False):
            outfile = os.path.join(outdir, label, name+".png")
            df = group.copy().sort_values("pos").reset_index(drop=True)
            fig = gen_pileup_fig(df, dataset_label=label, bgc=name.replace("_","."))
            
            #get annotation file:
            genbank_fp = os.path.join(genbank_dir, f"{name}.gbk")
            record = next(SeqIO.parse(genbank_fp, "genbank")) #Assumes only one record
            fig_anno = add_genomic_annotation(fig, record)
            fig.write_image(outfile, height=700, width=1000)
