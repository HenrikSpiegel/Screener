from pathlib import Path
import pandas as pd
import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--fuzzy_path", required=True, type=str, help="Fuzzy path for .counted files.")
parser.add_argument("-o", required=True, type=Path, help="Path for output directory")
args = parser.parse_args()

outdir = Path(args.o)
outdir.mkdir(parents=True, exist_ok=True)

path_glob = Path().glob(args.fuzzy_path)
grouped_dict = dict()
for path in path_glob:
    cat_name = path.stem
    if cat_name in grouped_dict:
        grouped_dict[cat_name].append(path)
    else:
        grouped_dict[cat_name] = [path]

regex_sample_name = re.compile(r"\/(\d+_?\d+GB)\/(sample_\d)\/")
for cat_name, cat_count_files in grouped_dict.items():
    dfs = []
    for count_file in cat_count_files:  
        regex_res = regex_sample_name.search(count_file.as_posix())
        sample_name = ".".join((regex_res.group(1),regex_res.group(2)))
        dfs.append(
            pd.read_csv(count_file, 
                     sep=" ", index_col=0, names=[sample_name])
        )
    df_cat = pd.concat(dfs, axis=1)
    df_cat.to_csv(outdir/(cat_name+".csv"))