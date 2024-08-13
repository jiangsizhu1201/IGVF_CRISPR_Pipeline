#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
import muon as mu
from mudata import MuData
import matplotlib.pyplot as plt
import os

def volcano(mdata: MuData):
    test_results = pd.DataFrame({k: v for k, v in mdata.uns['test_results'].items()})
    
    plt.scatter(x=test_results['log2_fc'], y=test_results['p_value'].apply(lambda x: -np.log10(x)), s=5, color="green", label="Not significant")
    
    # Highlight down or up-regulated genes
    down = test_results[(test_results['log2_fc'] <= -2) & (test_results['p_value'] <= 0.01)]
    up = test_results[(test_results['log2_fc'] >= 2) & (test_results['p_value'] <= 0.01)]
    
    plt.scatter(x=down['log2_fc'], y=down['p_value'].apply(lambda x: -np.log10(x)), s=10, label="Down-regulated", color="blue")
    plt.scatter(x=up['log2_fc'], y=up['p_value'].apply(lambda x: -np.log10(x)), s=10, label="Up-regulated", color="red")
    
    for index, row in down.iterrows():
        plt.text(x=row['log2_fc'], y=-np.log10(row['p_value']), s=row['gene_id'])

    low, high = plt.xlim()
    bound = max(abs(low), abs(high)) + 0.5
    
    plt.xlabel("log2 fold change")
    plt.ylabel("p value (log10)")
    plt.axvline(-2, color="grey", linestyle="--")
    plt.axvline(2, color="grey", linestyle="--")
    plt.axhline(2, color="grey", linestyle="--")
    plt.xlim(-bound, bound)
    plt.legend(loc="upper right")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a volcano plot from MuData")
    parser.add_argument("mdata_path", type=str, help="Path to the MuData file")
    
    args = parser.parse_args()
    
    # Load MuData
    mdata = mu.read(args.mdata_path)
    
    # Generate the volcano plot
    volcano(mdata)
    output_file = os.path.join("evaluation_output", "volcano_plot.png")

    plt.tight_layout()
    plt.savefig(output_file)
    plt.close() 
