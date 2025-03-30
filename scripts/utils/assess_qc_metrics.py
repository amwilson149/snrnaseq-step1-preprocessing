import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import sys
from matplotlib import pyplot as plt

# Helper functions to compute and visualize expression data
# metrics for QC, including n total molecules, n unique genes,
# and fraction of mitochondrial reads per barcode.

# The function to compute QC metrics  adds these measures to the
# input anndata object inplace.
def compute_qc_metrics(adata: ad.AnnData):
    # Calculate a few metrics and look at them
    if 'mito_frac' not in adata.obs.keys():
        sys.stderr.write('Calculating mito_frac, assuming mito genes start with\ngenome prefix plus "mt-" or "MT-"')
        # This works for our "<genome>________________MT-" format for mitochondrial genes in the
        # human case. Check that it also works for mouse genome data.
        if adata.var_names.str.contains('mt-',case=False).sum() == 0:
            sys.stderr.write('WARNING: no genes found that start with "mt-"')
        adata.obs['mito_frac'] = np.array(
            adata.X[:, adata.var_names.str.contains('mt-',case=False)].sum(axis=1)
        ).squeeze() / np.array(adata.X.sum(axis=1) + 1e-10).squeeze()
    adata.obs['n_umi'] = np.array(adata.X.sum(axis=1)).squeeze()
    adata.obs['n_gene'] = np.array((adata.X > 0).sum(axis=1)).squeeze()

def plot_qc_metric_distributions(adata: ad.AnnData,
        output_dir,
        metrics_to_plot=[
            'mito_frac',
            'n_umi',
            'n_gene'
            ]
        ):
    # Plot each of the metrics of interest
    for metric in metrics_to_plot:
        if metric in adata.obs.columns.values.tolist():
            # 04b.ii. Compute the upper percentiles of the
            # distribution of this metric and add this
            # information to a violin plot of the
            # metric's distribution
            plt.figure()
            ax = plt.subplot(1,1,1)
            sc.pl.violin(adata,
                    keys=metric,
                    show=False,
                    density_norm='width', # adjusted to accommodate new behavior
                    ax=ax
                    )
            x_tick_locs,x_tick_labels = plt.xticks()
            percentiles_to_plot = [90,95]
            y_pctiles = np.percentile(
                    adata.obs[metric].values.tolist(),
                    percentiles_to_plot)
            plt.plot([x_tick_locs[0]]*len(y_pctiles),
                    y_pctiles,
                    'x',
                    markersize=4,
                    linestyle='None')
            pctile_text = '_'.join([f'_'
                for _ in percentiles_to_plot
                ])
            plt.tight_layout()
            plt.savefig(f'{output_dir}/pre_filtering_pre_qc_{metric}_dist_all_w_pctiles_{pctile_text}.png',
                    dpi=300)
            plt.close()
        else:
            print(f'Metric {metric} not found in adata.obs. Skipping.')


