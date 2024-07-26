#!/usr/bin/env python

import argparse
import muon as mu
import anndata as ad
import pandas as pd
import numpy as np
from muon import MuData

import scanpy as sc
import scrublet as scr
import matplotlib.pyplot as plt
import os

def main():
    parser = argparse.ArgumentParser(description="Scrublet doublet detection for MuData object")
    parser.add_argument('--input', type=str, required=True, help='Path to the input MuData file')

    args = parser.parse_args()
    mdata = mu.read(args.input)

    adata = mdata.mod['gene']

    # Run scrublet
    scrub = scr.Scrublet(adata.X)
    adata.obs['doublet_scores'], adata.obs['predicted_doublets'] = scrub.scrub_doublets()
    scrub.plot_histogram()

    # Save plot
    if not os.path.exists('figures'):
        os.makedirs('figures')
        print(f"Directory '{'figures'}' created.")
    else:
        print(f"Directory already exists.")
    plt.savefig(f"figures/doublets_batch.png")
    
    print("Number of predicted doublets:", sum(adata.obs['predicted_doublets']))

    adata.obs['doublet_info'] = adata.obs["predicted_doublets"].astype(str)

    # Remove doublets
    adata = adata[adata.obs['predicted_doublets'] == False]

    # Save
    mdata.write('mdata_doublets.h5mu')

if __name__ == '__main__':
    main()
