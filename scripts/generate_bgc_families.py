
# Use an identity based clustering method going forward.

all_bgcs = ['NC_014328.1.region001', 'NC_014328.1.region002',
       'NC_014328.1.region003', 'NC_014328.1.region004',
       'NC_014328.1.region005', 'NZ_CP020566.1.region001',
       'NZ_CP053893.1.region001', 'NZ_CP053893.1.region002',
       'NZ_CP053893.1.region003', 'NZ_CP053893.1.region004',
       'NZ_CP053893.1.region005', 'NZ_CP053893.1.region006',
       'NZ_LT906445.1.region001', 'NZ_LT906470.1.region001',
       'NZ_LT906470.1.region002', 'NZ_LT906470.1.region003']

pre_assigned = ["NZ_LT906470.1.region003", "NZ_CP020566.1.region001", "NZ_LT906445.1.region001"]
manualclustermap = dict(
    ranthipeptide_group = ['NZ_CP020566.1.region001','NZ_LT906445.1.region001'],
    ranthipeptide_alone = ['NZ_LT906470.1.region003']
)
#all others are single member families.
manualclustermap.update({x:[x] for x in all_bgcs if x not in pre_assigned})

if __name__ == "__main__":
    import argparse
    import json
    from pathlib import Path
    import sys

    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--outfile", required=True, type=Path, help="location for json dump mapping families to member bgcs.")

    args = parser.parse_args()

    outfile = Path(args.outfile)
    outfile.parent.mkdir(parents=True, exist_ok=True)
    outfile.write_text(json.dumps(manualclustermap))

    print(f"Wrote cluster families to -> {outfile}", file=sys.stderr)