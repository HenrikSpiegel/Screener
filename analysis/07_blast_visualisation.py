from pathlib import Path
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly
from typing import List
from Bio import SeqIO

def get_sequence_lengths(mfasta):
    seq_len = {}
    for record in SeqIO.parse(mfasta, "fasta"):
        seq_len[record.name] = len(record.seq)
    return seq_len

def generate_underlying_scatter(query, subject) -> List[go.Scatter]:
    query_pos = "Q: "+query
    subject_pos = "S: "+subject
    base_traces = [
       go.Scatter({
         'line': {'color': 'lightgrey', 'dash': 'solid', 'width': 10},
         'mode': 'lines',
         'name': query,
         'orientation': 'h',
         'showlegend': True,
         'x': [    0, seq_lens[query]],
         'xaxis': 'x',
         'y': [query_pos, query_pos],
         'yaxis': 'y'
     }),
        go.Scatter({
         'line': {'color': 'lightgrey', 'dash': 'solid', 'width': 10},
         'mode': 'lines',
         'name': subject,
         'orientation': 'h',
         'showlegend': True,
         'x': [    0, seq_lens[subject]],
         'xaxis': 'x',
         'y': [subject_pos, subject_pos],
         'yaxis': 'y'
     }) 
    ]
    return base_traces

def generate_hsp_scatters(group) -> List[go.Scatter]:
    df = group.sort_values("qstart").reset_index(drop=True)
    scatters = []
    q = "Q: " + df.at[0, "qaccver"]
    s = "S: " + df.at[0, "saccver"]
    for row in df.itertuples():
        scatters.extend([
            go.Scatter(
                x=[row.qstart, row.qend],
                y=[q, q],
                line = {'color': plotly.colors.qualitative.Alphabet[row.Index], 'dash': 'solid', "width":5}
            ),
            go.Scatter(
                x=[row.sstart, row.send],
                y=[s, s],
                line = {'color': plotly.colors.qualitative.Alphabet[row.Index], 'dash': 'solid', "width":5}
            )
        ])
    return scatters


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", type=Path, required=True, help="fastafile containing input to blast")
    parser.add_argument("--blast", type=Path, required=True, help="blast output in outfmt 6")
    parser.add_argument("--similarity-table", type=Path, required=True, help="Table containing pairwise similarities for input seqs.")
    parser.add_argument("-o", type=Path, help="Overwrite default outout")
    parser.add_argument("-ph", type=int, default=2000, help="plot height")
    parser.add_argument("-pw", type=int, default=500, help="plot width")
    args = parser.parse_args()

    outdir = args.o or Path("results")/Path(__file__).stem
    outdir.mkdir(parents=True, exist_ok=True)
    sim_req = 5

    seq_lens = get_sequence_lengths(args.fasta)
    
    df_blast = pd.read_csv(args.blast, sep="\t")
    df_similarity = pd.read_csv(args.similarity_table, index_col=0,sep="\t")
    df_similarity.fillna(0, inplace=True)
    

    grouped = df_blast.query("qaccver != saccver").groupby(["qaccver","saccver"])
    grouped_similar = {name:grouped.get_group(name) for name in grouped.groups.keys() if df_similarity.loc[name]>5}

    subplot_titles = [
                        f"{name[0]} vs {name[1]}: Similarity {df_similarity.loc[name]:.0f}%" 
                        for name in grouped_similar             
                    ]
    fig = make_subplots(rows=len(grouped_similar), cols=1,
                    shared_xaxes=True,
                    subplot_titles=subplot_titles
                    )#,
                        #vertical_spacing=0.02)

    for i, (name, group) in enumerate(grouped_similar.items()):
        
        base_traces = generate_underlying_scatter(*name)
        fig.add_trace(base_traces[0], row=i+1, col=1)
        fig.add_trace(base_traces[1], row=i+1, col=1)
        for scat in generate_hsp_scatters(group):
            fig.add_trace(scat, row=i+1, col=1)
    fig.update_layout(height=2000, showlegend=False)
    fig.update_yaxes(title="")
    fig.update_layout(title="Visualization of BLAST results<br><sup> High Scoring Pairs (HSPs) shown as colored strecthes")

    fig.write_image(outdir/"Blast_comparison.png", height=args.ph, width=args.pw)









