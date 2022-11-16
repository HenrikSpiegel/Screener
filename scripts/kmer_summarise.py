import sys
import pandas as pd
import glob
import os
from typing import Sequence
import numpy as np
import statsmodels.api as sm

import configparser



def error_correction(values: Sequence[float], k:int=21, error_rate:float=0.03):
    correction_value = (1-error_rate)**k
    corrected_value = [x/correction_value for x in values]
    return corrected_value


def edgeloss_correction(values: Sequence[float], k:int=21, L:int=150):
    correction_value = 1 - ((k-1)/L)
    corrected_value = [x/correction_value for x in values]
    return corrected_value

def get_summary_frame(directory: str, file_fuzzy: str = "*.counted"):

    dfs = []
    for file in glob.glob(os.path.join(directory, file_fuzzy)):
        contig = file.rsplit("/",1)[1].split(".counted")[0]

        df = pd.read_csv(file, sep=" ", names=["mer","mer_count"])
        counts = df.loc[:, "count"]
        X = np.ones_like(counts)
        res = sm.NegativeBinomial(counts, X).fit()
        NB_mu = np.exp(res.params.const)
        
        #df["mer_count"] = 2*df["mer_count"]
        df_res = df.describe().transpose().reset_index(drop=True)
        df_res["Contig"] = contig
        df_res["count_nonzero"] = len(df.query("mer_count > 0"))
        df_res["count_unique"] = len(df.mer.drop_duplicates())
        df_res["count_unique_nonzero"] = len(df.query("mer_count > 0").mer.drop_duplicates())
        df_res["negBinom_mu"] = NB_mu
        df_res.rename(columns={"50%":"Depth median", "mean":"Depth avg"}, inplace=True)

        df_out = df_res[["Contig"]+[x for x in df_res.columns if x != "Contig"]]
        dfs.append(df_out)

    return pd.concat(dfs)



if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--directory",  required=True, help="Directory holding the count files from 'jellyfish query'.")
    parser.add_argument("--file_fuzzy", default="*.counted", help="Fuzzy filename for counted files.['*.counted']")

    parser.add_argument("-k", "--kmerlength", type=int, help="Overwrite default settings given in config/config")
    parser.add_argument("-l", "--readlength", type=float, help="Average Readlength of trimmed reads, Overwrite default settings given in config/config")
    parser.add_argument("-e", "--error-rate", type=float, help="Overwrite default settings given in config/config")

    parser.add_argument("-o", help="output file [--directory / kmer_summation.tsv]")

    args = parser.parse_args()

    if not args.o:
        args.o = os.path.join(args.directory, "kmer_summation.tsv")
    
    df = get_summary_frame(args.directory, args.file_fuzzy)
    if all([x == 0 for x in df["Depth median"]]):
        print("kmer_summarise.py: All 'Depth median' == 0", file=sys.stderr)
        #raise RuntimeError("kmer_summarise.py: All 'Depth median' == 0")

    config = configparser.ConfigParser()
    config.read("config/project_config.ini")


    kmerlength  = args.kmerlength   or config.getint("KmerQuantification", "KmerLength")
    read_length = args.readlength   or config.getint("KmerQuantification", "AverageReadLength")
    error_rate  = args.error_rate   or config.getfloat("KmerQuantification", "PerBaseErrorRate")

    df["Depth Error Corrected"]         = error_correction(df["Depth median"], k=kmerlength, error_rate=error_rate)
    df["Depth Error/Edge Corrected"]    = edgeloss_correction(df["Depth Error Corrected"], k=kmerlength, L=read_length)

    df.to_csv(args.o, index=False, sep="\t")

