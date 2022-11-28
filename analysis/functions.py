from typing import List
import numpy as np

def add_genomic_annotation(fig, genbank_record, custom_annot: List[dict]=[]):
    y1_relative = 0.12
    anno_yshift_relative = 0.85
    color_map = {
        "CDS":"grey",
        "rRNA": "darkred"
    }

    y_magnitude = fig.data[0]["y"].max()
    shape_y1 = -(y1_relative * y_magnitude)
    anno_y1  = -(y1_relative * y_magnitude)
    anno_yshift = -(anno_yshift_relative*y_magnitude)

    for feat in genbank_record.features:
        if feat.type=="region":
            pass
        if feat.type=="CDS":
            location = (min(feat.location), max(feat.location))
            locus_tag = feat.qualifiers["locus_tag"]
            fig.add_shape(
                name=locus_tag[0],
                x0 = location[0],
                x1 = location[1],
                y0=-0.1,
                y1=shape_y1,
                yref="y",
                fillcolor=color_map[feat.type]
            )
            fig.add_annotation(
                text = locus_tag[0],
                x = np.mean(location),
                y=anno_y1,
                yshift=-35,
                showarrow=False,
                textangle=90,
                align="right"
            )
    for entry in custom_annot: #does nothing if custom_annot == []
        fig.add_shape(
                name=entry["locus_tag"],
                x0 = entry["location"][0],
                x1 = entry["location"][1],
                y0=-0.1,
                y1=shape_y1,
                yref="y",
                fillcolor=color_map[entry["type"]]
        )
        fig.add_annotation(
            text = entry["locus_tag"],
            x = np.mean(entry["location"]),
            y=anno_y1,
            yshift=-35,
            showarrow=False,
            textangle=90,
            align="right"
        )
    return fig