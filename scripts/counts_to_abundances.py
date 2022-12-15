from pathlib import Path
import json
import pandas as pd
import statsmodels.api as sm
import numpy as np

def error_correction(values, k:int=21, error_rate:float=0.03):
    correction_value = (1-error_rate)**k
    corrected_value = [x/correction_value for x in values]
    return corrected_value


def edgeloss_correction(values, k:int=21, L:int=150):
    correction_value = 1 - ((k-1)/L)
    corrected_value = [x/correction_value for x in values]
    return corrected_value

def combined_correction(values, k, L, error_rate):
    err_corrected = error_correction(values=values, k=k, error_rate=error_rate)

    return edgeloss_correction(
        values=err_corrected,
        k=k, L=L
    )

def negbinom_mu(column:pd.Series):
    nb_fit = sm.NegativeBinomial(column, np.ones_like(column)).fit(disp=0, start_params=[1,1]) #disp=0 == quiet
    nb_param_mu = np.exp(nb_fit.params.const)
    return nb_param_mu

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--count-matrices", required=True, type=Path, help="dir containing count matrices")
    parser.add_argument("--mag-screened", required=True, type=Path, help="Dir containing extracted flat _kmers.tsv files")
    parser.add_argument("--read-len-json", required=True, type=Path, help="file containing avg readlenghts for the samples.")
    parser.add_argument("--kmer-len", required=True, type=int)
    parser.add_argument("--error-rate-est", required=True, type=float)
    parser.add_argument("-o", required=True, type=Path)
    args = parser.parse_args()

    ###

    dir_count_matrices = args.count_matrices       #WD_DATA / "kmer_quantification/count_matrices"
    dir_MAGinator      = args.mag_screened                         #WD_DATA / "MAGinator/screened_flat"

    outdir = args.o
    outdir.mkdir(parents=True, exist_ok=True)

    json_avg_readlen = args.read_len_json
    error_rate_estimate = args.error_rate_est
    kmer_len = args.kmer_len

    ####

    dict_avg_readlen = json.loads(json_avg_readlen.read_text())
    sample_to_readlen = {
        ".".join((Path(k).parent.name, Path(k).name)):v
        for k,v in dict_avg_readlen.items()
    }

    count_matrices = [x for x in dir_count_matrices.glob("*.tsv") if not x.name.startswith("counts_all")]

    if not count_matrices:
        raise FileNotFoundError(f"dir_count_matrices/*.tsv")

    abundances = []
    for fp_count in count_matrices:
        catalogue_name = fp_count.stem
        fp_catalogue = dir_MAGinator/ (catalogue_name+"_kmers.csv")
        df_catalogue_MAG = pd.read_csv(fp_catalogue).melt(var_name = "method", value_name="kmer")
        
        
        df_count     = pd.read_csv(fp_count, index_col=0, sep="\t")
        df_count.reset_index(inplace=True)
        df_count.rename(columns = {'index':'kmer'}, inplace=True)
        
        #add the raw input to catalouge:
        df_catalogue = pd.concat([
            df_catalogue_MAG,
            pd.DataFrame({'method': f'raw_{len(df_count)}', 'kmer':df_count.kmer})
        ])
        
        # merge with counts
        df_combined = df_catalogue.merge(df_count, on="kmer", how="left")
        df_combined_long = df_combined.melt(id_vars=["method","kmer"], var_name="dataset_sample", value_name="count")
        
        # aggregate to estimates
        df_estimates = df_combined_long\
            .groupby(["method","dataset_sample"])["count"]\
            .agg([negbinom_mu, np.median])\
            .reset_index()
        
        # Apply corrections
        corrected = []
        for ds_sample, df_i in df_estimates.groupby("dataset_sample"):
            avg_readlen = sample_to_readlen[ds_sample]

            for est_col in ["negbinom_mu","median"]:
                df_i[est_col+"_corr"] = combined_correction(
                    df_i[est_col].values, 
                    k=kmer_len, L=avg_readlen, error_rate=error_rate_estimate)
            corrected.append(df_i)
        df_corrected = pd.concat(corrected).reset_index(drop=True)
        
        df_corrected["catalogue_name"] = catalogue_name
        
        outfile = outdir / (catalogue_name+"_abundances.csv")
        df_corrected.to_csv(outfile, index=False)
        
        abundances.append(df_corrected)

    df_abundances = pd.concat(abundances) 
    outfile = outdir / ("total_abundances.csv")
    df_abundances.to_csv(outfile, index=False)