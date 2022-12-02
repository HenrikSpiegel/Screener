import numpy as np
from pathlib import Path
from Bio import SeqIO

def prep_blast_demo(seq_file, outdir, n_shuffled_sequences=2, n_chunks=5) -> Path:
    data_dir = Path(outdir)
    data_dir.mkdir(parents=True, exist_ok=True)
    record = next(SeqIO.parse(seq_file, "fasta"))
    sequence = record.seq.__str__()
    name = "SequenceX"

    chunk_size = int((len(sequence)/n_chunks) + 0.5) #ceil
    sequence_chunks = [sequence[i*chunk_size:(i+1)*chunk_size] for i in range(n_chunks)]
    shuffled_seqs = ["".join(np.random.choice(sequence_chunks, size=n_chunks, replace=False)) for i in range(n_shuffled_sequences)]

    outfile = data_dir / 'demo.fa'
    entries = [f">{name}\n{sequence}"]+[f">{name}-({i+1})\n{seq}" for i, seq in enumerate(shuffled_seqs)]
    outfile.write_text("\n".join(entries))
    return data_dir, outfile

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--n_shuffled_sequences", type=int)
    parser.add_argument("--n_chunks", type=int)
    parser.add_argument("--seq-file", type=Path, required=True)
    parser.add_argument("-o", type=Path, required=True)
    args = parser.parse_args()

    prep_blast_demo(
        seq_file=args.seq_file, 
        outdir=args.o, 
        n_shuffled_sequences=args.n_shuffled_sequences, 
        n_chunks=args.n_chunks)