from pathlib import Path
import json
import sys

def mcl2json(mcl_file):
    clusters_str = mcl_file.read_text().strip().split("\n")

    cluster_dict = {
        f"cluster_{i:02d}":[member for member in cluster_line.split()] 
        for i, cluster_line in enumerate(clusters_str)
    }
    return json.dumps(cluster_dict,indent=True)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--mcl-files", type=Path, nargs="+", help="list of MCL files")
    parser.add_argument("-o", type=Path, help="direct output files to another dir.")
    args = parser.parse_args()

    print("parsing:\n"+"\n".join(x.as_posix() for x in args.mcl_files), file=sys.stderr)

    for file in args.mcl_files:
        mcl_json = mcl2json(file)

        outname = file.name+".json"
        outdir = args.o or file.parent

        fp_out = outdir / outname
        fp_out.write_text(mcl_json)