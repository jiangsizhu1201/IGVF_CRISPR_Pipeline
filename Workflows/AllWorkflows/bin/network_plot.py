#!/usr/bin/env python
import os
import argparse
import pandas as pd
import numpy as np
import muon as mu
from mudata import MuData

import networkx as nx
import igraph 

import matplotlib.pyplot as plt

from typing import Optional

def add_nontargeting_gene_exp_col(
    mdata: MuData,
):
    # Get the guides that were non-targeting
    nontargeting_mask = mdata.mod["guide"].var["intended_target_name"] == "non-targeting"
    targeting_mask = ~nontargeting_mask

    # Get the counts across barcodes for the non-targeting guides and the targeting guides
    nontargeting_counts = mdata.mod["guide"][:, nontargeting_mask].copy()
    targeting_counts = mdata.mod["guide"][:, targeting_mask].copy()

    # Determine the barcodes that recieved at least one non-targeting guides (non-zero counts) AND no targeting guides (zero counts)
    cells_with_nontargeting_mask = nontargeting_counts.obs.index[(np.sum(nontargeting_counts.X, axis=1) > 0).A.ravel()]
    cells_with_no_targeting_mask = targeting_counts.obs.index[(np.sum(targeting_counts.X, axis=1) == 0).A.ravel()]
    cells_with_only_nontargeting = np.intersect1d(cells_with_nontargeting_mask, cells_with_no_targeting_mask)
    
    # If no cells have only non-targeting guides, handle the situation
    if len(cells_with_only_nontargeting) == 0:
        print("Warning: No cells with only non-targeting guides found. Skipping mean calculation.")
        mdata.mod["gene"].var["mean_nontargeting_expression"] = np.nan
        return
    
    # 
    nontargeting_adata = mdata.mod["gene"][cells_with_only_nontargeting].copy()
    
    # Get gene expression pd.Series where the index is the gene name and the values are the mean expression across cells
    gene_expression = nontargeting_adata.X.mean(axis=0).A1
    mdata.mod["gene"].var["mean_nontargeting_expression"] = gene_expression

def plot_network(
    mdata: MuData,
    central_node: Optional[str] = None,
    source_column: str = "source",
    target_column: str = "target",
    weight_column: Optional[str] = None,
    min_weight: Optional[float] = None,
    node_size_column: Optional[str] = None,
    results_key: Optional[str] = "test_results",
):
    # Get and filter the results dataframe
    results_df = pd.DataFrame({k: v for k, v in mdata.uns[results_key].items()})
    if min_weight is not None:
        results_df = results_df[results_df[weight_column].abs() >= min_weight]
    if central_node is not None:
        results_df = results_df[results_df[source_column] == central_node]

    # Create the network 
    G = nx.DiGraph()
    for i, row in results_df.iterrows():
        G.add_edge(row[source_column], row[target_column], weight=row[weight_column])
    pos = nx.circular_layout(G)

    # Draw nodes based on node size column if provided
    if node_size_column is not None:
        node_size = results_df.set_index(target_column)[node_size_column].to_dict()
        sizes = [node_size.get(n, 10) for n in G.nodes]
        sizes = [s * 100 for s in sizes]
        nx.draw(G, pos, node_color="skyblue", node_size=sizes, edge_cmap=plt.cm.Blues, arrowsize=20)
    else:
        nx.draw(G, pos, with_labels=True, node_color="skyblue", node_size=150, edge_cmap=plt.cm.Blues, arrowsize=20)

    # Offset labels
    if central_node is not None:
        label_pos = {k: (v[0]+0.2, v[1] + 0.05) for k, v in pos.items() if k != central_node}  # Adjust 0.1 as needed for the offset
        label_pos[central_node] = pos[central_node]
        nx.draw_networkx_labels(G, label_pos, font_size=10)
    else:
        label_pos = {k: (v[0]+0.2, v[1] + 0.05) for k, v in pos.items()}
        nx.draw_networkx_labels(G, pos, font_size=8)
    
    # Draw edges
    edge_labels = {(u, v): f"{d['weight']:.2f}" for u, v, d in G.edges(data=True)}
    nx.draw_networkx_edges(G, pos, edge_color=[d["weight"] for u, v, d in G.edges(data=True)], edge_cmap=plt.cm.coolwarm, arrowsize=5)
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels)
    
    # Plt
    output_dir = "evaluation_output"
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, "network_plot.png")

    plt.tight_layout()
    plt.savefig(output_file)
    plt.close() 

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot a network from MuData")
    parser.add_argument("mdata_path", type=str, help="Path to the MuData file")
    parser.add_argument("--central_node", type=str, default=None, help="Central node to plot (optional)")
    parser.add_argument("--source_column", type=str, default="intended_target_name", help="Source column name in test results")
    parser.add_argument("--target_column", type=str, default="gene_id", help="Target column name in test results")
    parser.add_argument("--weight_column", type=str, default="log2_fc", help="Weight column name in test results (optional)")
    parser.add_argument("--min_weight", type=float, default=None, help="Minimum weight to filter edges (optional)")
    parser.add_argument("--node_size_column", type=str, default="p_value", help="Column to determine node size (optional)")
    parser.add_argument("--results_key", type=str, default="test_results", help="Key for test results in mdata.uns")
    
    args = parser.parse_args()
    
    # Load MuData
    mdata = mu.read(args.mdata_path)
    
    # Plot network
    plot_network(
        mdata,
        central_node=args.central_node,
        source_column=args.source_column,
        target_column=args.target_column,
        weight_column=args.weight_column,
        min_weight=args.min_weight,
        node_size_column=args.node_size_column,
        results_key=args.results_key,
    )
