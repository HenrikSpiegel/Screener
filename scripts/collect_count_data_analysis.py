from pathlib import Path
import pandas as pd

# Collect initial 10k counts.
wd = Path("data/ibdmdb")
dir_kmer_count_matrices = wd / "kmer_quantification/count_matrices"
count_fuzzy = wd.glob("kmer_quantification/SR*/counts/*.counted")


grouped_dict = dict()
for path in wd.glob("kmer_quantification/SR*/counts/*.counted"):
    cat_name = path.stem
    if cat_name in grouped_dict:
        grouped_dict[cat_name].append(path)
    else:
        grouped_dict[cat_name] = [path]

all_counts = []
for cluster_name, cluster_count_files in grouped_dict.items():
    df_cluster_count = pd.concat(
        (pd.read_csv(fp_count, index_col=0, names=[fp_count.parent.parent.name], sep=" ")
        for fp_count in cluster_count_files),
        axis=1
    )
    df_cluster_count.to_csv(dir_kmer_count_matrices/(cluster_name+".tsv"), sep="\t")
    all_counts.append(df_cluster_count)

pd.concat(all_counts).to_csv(dir_kmer_count_matrices/"counts_all.tsv", sep="\t")