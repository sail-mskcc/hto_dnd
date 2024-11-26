import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

def compare_distributions(raw_data, dsb_data, output_path):
    n_proteins = raw_data.shape[1]
    fig, axes = plt.subplots(n_proteins, 2, figsize=(15, 4*n_proteins))
    
    for i, protein in enumerate(raw_data.columns):
        # Raw data
        raw_data_protein = raw_data[protein]
        dsb_data_protein = dsb_data[protein]
        
        # Calculate x-axis limits for raw data
        raw_xlim = np.percentile(raw_data_protein, [0, 99])
        
        sns.histplot(raw_data_protein, ax=axes[i, 0], kde=True)
        axes[i, 0].set_title(f'{protein} - Raw')
        axes[i, 0].set_xlim(raw_xlim)
        
        # DSB data
        # Calculate x-axis limits for DSB data
        dsb_xlim = np.percentile(dsb_data_protein, [0, 99.5])
        
        sns.histplot(dsb_data_protein, ax=axes[i, 1], kde=True)
        axes[i, 1].set_title(f'{protein} - DSB')
        axes[i, 1].set_xlim(dsb_xlim)
        
        # Add text with min, max, and median values
        axes[i, 0].text(0.05, 0.95, f"Min: {raw_data_protein.min():.2f}\nMax: {raw_data_protein.max():.2f}\nMedian: {raw_data_protein.median():.2f}", 
                        transform=axes[i, 0].transAxes, verticalalignment='top')
        axes[i, 1].text(0.05, 0.95, f"Min: {dsb_data_protein.min():.2f}\nMax: {dsb_data_protein.max():.2f}\nMedian: {dsb_data_protein.median():.2f}", 
                        transform=axes[i, 1].transAxes, verticalalignment='top')
    
    plt.tight_layout()
    
    # Save the plot
    plt.savefig(output_path)
    plt.close()

def create_visualization(adata, output_path):
    # Extract raw and DSB-normalized data
    raw_data = pd.DataFrame(adata.X, columns=adata.var_names)
    dsb_data = pd.DataFrame(adata.layers['dsb_normalized'], columns=adata.var_names)
    
    # Create the output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Create and save the plot
    compare_distributions(raw_data, dsb_data, output_path)