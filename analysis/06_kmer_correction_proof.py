import random
import sys
import statsmodels.api as sm
import pandas as pd
from tqdm import trange
import plotly.express as px
import numpy as np
from pathlib import Path

import plotly.io as pio
pio.kaleido.scope.mathjax = None
#
output_dir = Path("results/") / Path(__file__).stem

# simulation controls:
experiments = 10000
error_rates = np.array([0,1,2,3,4,5])
read_lenghts = np.array([50,100,150,200,300])
k = 21

catalogue_size = 20

def kmerize(seq, k):
    return np.array(["".join(seq[i:i+k]) for i in range(len(seq)-k+1)])

bases = "ATCG"
bases_set={x for x in bases}

substituion_dict = {b:[*(bases_set).difference(set(b))] for  b in bases_set}
def substitute(base):
    return random.choice(substituion_dict[base])

def error_rate_variations():
    # run simmulations
    print("Running error rate simulations", file=sys.stderr)

    L = 150
    sequence_L = 3*L

    count = np.zeros((len(error_rates), experiments))
    for j in trange(len(error_rates)):
        err_rate = error_rates[j]
        for i in trange(experiments, leave=False
                    ):
            true_sequence = np.array([random.choice(bases) for i in range(sequence_L)])
            true_read     = true_sequence[0:L]
            catalogue = set(random.choice(kmerize(true_sequence, k)) for i in range(catalogue_size))


            read = true_read.copy()
            is_error = np.random.randint(1,100, L) <= err_rate
            read[is_error] = [substitute(x) for x in read[is_error]]

            kmers_obs = list(kmerize(read, k))
            count[j, i] = sum(x in catalogue for x in kmers_obs) / len(catalogue) #TODO: Use a distribution based center metric (NB center)

    print("Adding error corrections and plotting", file=sys.stderr)
    # Add error correction
    error_rate_correction = (1-error_rates/100)**k
    count_error_corrected = count / error_rate_correction[:,np.newaxis]

    edge_loss_correction = 1 - ((k-1)/L)
    count_fully_corrected = count_error_corrected / edge_loss_correction

    cols = [f"{x}%" for x in error_rates]
    df_count = pd.DataFrame(count.T, columns = cols)
    df_count["type"] = "count"

    df_count_err = pd.DataFrame(count_error_corrected.T, columns = cols)
    df_count_err["type"] = "Error corrected count"

    df_count_final = pd.DataFrame(count_fully_corrected.T, columns = cols)
    df_count_final["type"] = "Error and Edge loss corrected count"

    print(f"Writing results to -> {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    df_total = pd.concat([df_count, df_count_err, df_count_final])
    df_total.to_csv(output_dir/"counts_err.csv", index=False)

    df_total_long = df_total.melt(id_vars="type", var_name="Error Rate", value_name="Coverage Estimate")
    fig = px.box(df_total_long, x="type", y="Coverage Estimate", color="Error Rate",
            title="Coverage estimates from simulated read of various sequencing errors<br><sup>Simulated coverage is 33% with read length 150",
            labels={"type":""})
    fig.add_hline(y=1/3, line_dash='dot', annotation_text="Coverage", annotation_position="top left")
    fig.update_layout(template="plotly_white")
    fig.write_image(output_dir/"coverage_estimates_error.png", width=1000, height=500)

def read_length_variations():
    # run simmulations
    print("Running read length simulations", file=sys.stderr)
    err_rate = 3
    count = np.zeros((len(read_lenghts), experiments))
    for j in trange(len(read_lenghts)):
        L = read_lenghts[j]
        sequence_L = 3*L
        for i in trange(experiments, leave=False):
            true_sequence = np.array([random.choice(bases) for i in range(sequence_L)])
            true_read     = true_sequence[0:L]
            catalogue = set(random.choice(kmerize(true_sequence, k)) for i in range(catalogue_size))

            read = true_read.copy()
            is_error = np.random.randint(1,100, L) <= err_rate
            read[is_error] = [substitute(x) for x in read[is_error]]

            kmers_obs = list(kmerize(read, k))
            count[j, i] = sum(x in catalogue for x in kmers_obs) / len(catalogue) #TODO: Use a distribution based center metric (NB center)

    print("Adding error corrections and plotting", file=sys.stderr)
    # Add error correction
    error_rate_correction = (1-err_rate/100)**k
    count_error_corrected = count / error_rate_correction

    edge_loss_correction = 1 - ((k-1)/read_lenghts)
    count_fully_corrected = count_error_corrected / edge_loss_correction[:,np.newaxis]

    cols = [f"{x}bp" for x in read_lenghts]
    df_count = pd.DataFrame(count.T, columns = cols)
    df_count["type"] = "count"

    df_count_err = pd.DataFrame(count_error_corrected.T, columns = cols)
    df_count_err["type"] = "Error corrected count"

    df_count_final = pd.DataFrame(count_fully_corrected.T, columns = cols)
    df_count_final["type"] = "Error and Edge loss corrected count"

    print(f"Writing results to -> {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    df_total = pd.concat([df_count, df_count_err, df_count_final])
    df_total.to_csv(output_dir/"counts_readlength.csv", index=False)

    df_total_long = df_total.melt(id_vars="type", var_name="Read Length", value_name="Coverage Estimate")
    fig = px.box(df_total_long, x="type", y="Coverage Estimate", color="Read Length",
            title="Coverage estimates from simulated read of various read lengths<br><sup>Simulated coverage is 33% at 3% sequencing error",
            labels={"type":""})
    fig.add_hline(y=1/3, line_dash='dot', annotation_text="Coverage", annotation_position="top left")
    fig.update_layout(template="plotly_white")
    fig.write_image(output_dir/"coverage_estimates_readlength.png", width=1000, height=500)

if __name__ == "__main__":
    print("Running" + __file__, file=sys.stderr)
    error_rate_variations()
    read_length_variations()
    print("Finished: "+__file__, file=sys.stderr)