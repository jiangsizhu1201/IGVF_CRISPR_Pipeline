#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
from collections import defaultdict
import muon as mu
import os

def igv(mdata):
    coordinate_dict = {}

    for index, row in mdata.mod["gene"].var.iterrows():
        if np.isnan(row["gene_start"]) or np.isnan(row["gene_end"]):
            continue
        coordinate_dict[index] = [row["gene_chr"], row["gene_start"], row["gene_end"]]

    for index, row in mdata.mod["guide"].var.iterrows():
        if row["intended_target_name"] in coordinate_dict or row["intended_target_name"] == "non-targeting":
            continue
        coordinate_dict[row["intended_target_name"]] = [row["intended_target_chr"], row["intended_target_start"], row["intended_target_end"]]

    bedpe = defaultdict(list)
    bedgraph = defaultdict(list)
    test_results = pd.DataFrame({k: v for k, v in mdata.uns['test_results'].items()})

    for index, row in test_results.iterrows():
        if row["intended_target_name"] == row["gene_id"]:
            # PROMOTER
            bedgraph["chr"].append(coordinate_dict[row["intended_target_name"]][0])
            bedgraph["start"].append(coordinate_dict[row["intended_target_name"]][1])
            bedgraph["end"].append(coordinate_dict[row["intended_target_name"]][2])
            bedgraph["posterior_probability"].append(row["posterior_probability"])
            bedgraph["log2_fc"].append(row["log2_fc"])
        else:
            bedpe["chr1"].append(coordinate_dict[row["intended_target_name"]][0])
            bedpe["start1"].append(coordinate_dict[row["intended_target_name"]][1])
            bedpe["end1"].append(coordinate_dict[row["intended_target_name"]][2])

            bedpe["chr2"].append(coordinate_dict[row["gene_id"]][0])
            bedpe["start2"].append(coordinate_dict[row["gene_id"]][1])
            bedpe["end2"].append(coordinate_dict[row["gene_id"]][2])

            bedpe["p_value"].append(row["p_value"])
            bedpe["log2_fc"].append(row["log2_fc"])

    bedpe_df = pd.DataFrame(bedpe)
    bedgraph_df = pd.DataFrame(bedgraph)
    
    return bedpe_df, bedgraph_df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process MuData and generate bedpe and bedgraph files")
    parser.add_argument("mdata_path", type=str, help="Path to the MuData file")
    parser.add_argument("--bedpe_file", type=str, default="output.bedpe", help="Filename for the bedpe output (default: output.bedpe)")
    parser.add_argument("--bedgraph_file", type=str, default="output.bedgraph", help="Filename for the bedgraph output (default: output.bedgraph)")
    
    args = parser.parse_args()
    
    # Load MuData
    mdata = mu.read(args.mdata_path)
    
    # Process the data
    bedpe_df, bedgraph_df = igv(mdata)
    

    # Save the output files
    bedpe_path = os.path.join("evaluation_output", args.bedpe_file)
    bedgraph_path = os.path.join("evaluation_output", args.bedgraph_file)
    
    bedpe_df.to_csv(bedpe_path, sep="\t", index=False, header=False)
    bedgraph_df.to_csv(bedgraph_path, sep="\t", index=False, header=False)
    
    print(f"bedpe file saved to {bedpe_path}")
    print(f"bedgraph file saved to {bedgraph_path}")
