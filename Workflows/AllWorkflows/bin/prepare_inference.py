#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
import muon as mu

def main(guide_inference, mudata_path):
    # read in files
    guide_inference = pd.read_csv(guide_inference)
    mudata = mu.read(mudata_path)

    # adding 'pairs_to_test' to mudata uns
    ## remove the gene_id that not included in the mudata
    include = set(guide_inference['gene_name']).intersection(set(mudata.mod['gene'].var.index))
    subset = guide_inference[guide_inference['gene_name'].isin(include)]
    mudata.uns['pairs_to_test'] = subset.to_dict(orient='list')

    ## rename uns 
    key_mapping = {'gene_name': 'gene_id',
                'guide_id': 'guide_id',
                'intended_target_name': 'intended_target_name',
                'pair_type': 'pair_type'}

    mudata.uns['pairs_to_test'] = {key_mapping[k]: v for k, v in mudata.uns['pairs_to_test'].items()}
    
    ## save the mudata
    mudata.write("mudata_inference_input.h5mu")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process mudata and add guide_inference to uns.")
    parser.add_argument('guide_inference', type=str, help='Path to the input inference file')
    parser.add_argument('mudata_path', type=str, help='Path to the input MuData file')

    args = parser.parse_args()
    main(args.guide_inference, args.mudata_path)
