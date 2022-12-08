from pathlib import Path
from typing import List
import configparser
import argparse
import pandas as pd

def get_simulation_overview(simulation_overview_files:List[Path], camisim_config_file:Path):
    dfs = []
    for summary_file in simulation_overview_files:
        df = pd.read_csv(summary_file, sep=",")
        df["dataset"] = summary_file.parent.name
        dfs.append(df)
    df_simulation = pd.concat(dfs)
    
    #Grab config data:
    camisim_config = configparser.ConfigParser()
    #Here we assume only size (ie. readdepth) changes.
    config_file = Path(camisim_config_file)
    camisim_config.read(config_file)
    df_simulation["err_type"]=camisim_config.get("ReadSimulator","type")
    try:
        df_simulation["err_profile"]=camisim_config.getfloat("ReadSimulator","profile")
    except:
        df_simulation["err_profile"]=camisim_config.get("ReadSimulator","profile")
    df_simulation["fragment_size"]=camisim_config.getint("ReadSimulator","fragments_size_mean")
    
    return df_simulation

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--camisim-overview-files", nargs="+", type=Path, required=True)
    parser.add_argument("--camisim-config", type=Path, required=True)
    parser.add_argument("--outfile", type=Path, required=True)
    args = parser.parse_args()

    df = get_simulation_overview(simulation_overview_files=args.camisim_overview_files, camisim_config_file=args.camisim_config)
    df.to_csv(args.outfile, sep="\t", index=False)