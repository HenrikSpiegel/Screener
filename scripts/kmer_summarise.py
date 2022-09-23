import pandas as pd
import glob
import os

def get_summary_frame(directory: str, file_fuzzy: str = "*.counted"):

    dfs = []
    for file in glob.glob(os.path.join(directory, file_fuzzy)):
        contig = file.rsplit("/",1)[1].split(".counted")[0]

        df = pd.read_csv(file, sep=" ", names=["mer","mer_count"])
        df_res = df.describe().transpose().reset_index(drop=True)
        df_res["Contig"] = contig
        df_res["count_nonzero"] = len(df.query("mer_count > 0"))
        df_res["count_unique"] = len(df.mer.drop_duplicates())
        df_res["count_unique_nonzero"] = len(df.query("mer_count > 0").mer.drop_duplicates())
        df_res.rename(columns={"50%":"Depth median", "mean":"Depth avg"}, inplace=True)
        df_out = df_res[["Contig"]+[x for x in df_res.columns if x != "Contig"]]
        dfs.append(df_out)

    return pd.concat(dfs)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--directory",  required=True, help="Directory holding the count files from 'jellyfish query'.")
    parser.add_argument("--file_fuzzy", default="*.counted", help="Fuzzy filename for counted files.['*.counted']")
    parser.add_argument("-o", help="output file [--directory / kmer_summation.tsv]")

    args = parser.parse_args()

    if not args.o:
        args.o = os.path.join(args.directory, "kmer_summation.tsv")
    
    df = get_summary_frame(args.directory, args.file_fuzzy)
    df.to_csv(args.o, index=False, sep="\t")

