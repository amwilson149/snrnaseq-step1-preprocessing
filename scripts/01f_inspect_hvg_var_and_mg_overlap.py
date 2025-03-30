import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import os
import sys
import argparse
import yaml
from utils.consts import *

# This script was used to select the
# number of HVGs to retain during
# expression data preprocessing.

# It reads in expression data that
# has undergone only the heuristic
# QC filtering step of preprocessing
# but no HVG-based filtering (the pre-
# Harmony expression data)
# and it produces a plot of rank-ordered
# gene variability (that is, HVG
# rankings).

# This script also computes the fraction
# of cell type marker genes from our database
# present among the top n HVGs,
# as a function of n, to illustrate how the
# number of cell type markers increases as
# the number of retained HVGs increases,
# enabling a choice about how many HVGs
# can be retained to increase the number
# of cell type marker genes in the data
# but avoid including many low-variability
# genes in the expression data.

# This script also takes a list of genes that are indicators
# of inflammation and a list of marker genes, and
# uses the Harmony-corrected, completely processed
# data to plot the expression if they are ranked
# to be within the set of HVGs that were retained
# in downstream analysis.

# This script uses modified versions of scaling and other functions that
# work with data types smaller than int64 and float64, so that a larger
# number of HVGs can be retained in the analysis.

# Script setup

# 00. Create argparse object for reading input arguments
parser = argparse.ArgumentParser(description='QC and preprocessing for aggregated, raw expression data.')
parser.add_argument(f'--config-yaml-path',type=str,required=True)
args = parser.parse_args()

# 01. Get analysis parameters from configuration file
cfg = None
with open(args.config_yaml_path) as f:
    cfg = yaml.safe_load(f)
# 01a. Dataset information
# 01a.i. Data root directory
qc_root_dir = cfg.get('root_dir')
# 01a.ii. Data directory with preprocessed expression data files
preprocessed_data_dir = qc_root_dir + f'/' + cfg.get('preprocessed_data_dir')
# 01a.iii. Input expression data file name
input_data_fn = cfg.get('dbl_rm_agg_raw_expr_data_fn')
# 01a.iv. QCed output expression data file name
output_data_no_batch_corr_fn = cfg.get('preproc_exp_data_no_batch_corr_raw_fn')
# 01a.v. Marker gene database directory
mg_database_dir = qc_root_dir + f'/' + cfg.get('mg_db_dir')
# 01a.vi. Marker gene database file name
mg_database_fn = cfg.get('mg_db_fn')
# 01a.vii. Set up output directory
# for HVG and marker gene overlap plots
hvg_mg_out_dir = f'{preprocessed_data_dir}/hvg_ranked_var_mg_overlaps'
if not os.path.isdir(hvg_mg_out_dir):
    os.system(f'mkdir -p {hvg_mg_out_dir}')
# 01a.viii. QCed, batch-corrected expression data directory
finalized_expr_data_dir = qc_root_dir + f'/' + cfg.get('post_qc_preproc_data_dir')
# 01a.ix. QCed, batch-corrected expression data
# file name
fn_expr_data_fn = cfg.get('post_qc_preproc_data_fn')
# 01a.x. QC hyperparameter dictionary
qc_hyperparam_dict = cfg.get('qc_hyperparams')
N_HVGS = qc_hyperparam_dict['n_hvgs']

# 01b. Preprocessing parameters
# 01b.i. Gene name prefix (to remove)
gene_prefix = 'GRCh38_______________'

# 02. Import expression data and marker
# gene database
# 02a. Import expression data
ahc_pre_harm_fn = f'{preprocessed_data_dir}/{output_data_no_batch_corr_fn}'
ahc = ad.read_h5ad(ahc_pre_harm_fn)

# 02b. Import fully QCed and preprocessed expression data
# (but which still has not been cell typed)
ahc_ct_fn = f'{finalized_expr_data_dir}/{fn_expr_data_fn}'
ahc_ct = ad.read_h5ad(ahc_ct_fn)

# 02c. Import marker gene database
mg_db_full_fn = f'{mg_database_dir}/{mg_database_fn}'
mg_db = pd.read_csv(
        mg_db_full_fn,
        index_col=0
        )

# 03. Get the .var matrix for the pre-processed 
# expression data, which has gene-specific
# annotations including HVG variances and rankings
ahc_var_cp = ahc.var.copy()
# 03a. Delete the full expression data object
# to save space
del ahc

# 03b. Pull HVG information
# 03b.i. HVG rankings (this variable is only
# assigned up to the number of HVGs to keep that were
# specified)
hvg_rank = ahc_var_cp['highly_variable_rank'].values.tolist()
# 03b.ii. Normed variances used for ranking genes
# and selecting HVGs
hvg_var_norm = ahc_var_cp['variances_norm'].values.tolist()
# 03b.iii. Get the list of gene names that go along with
# these rankings
gene_names_raw = ahc_var_cp.index.values.tolist()
# 03b.iv. Get a reformatted version of the gene name
# list with the prefix removed
gene_names = [
        _.split(gene_prefix)[-1]
        for _ in
        gene_names_raw
        ]

# 04. Sort the top HVGs according to their actual rankings
# (which preserves the actual tie-breaking done in our data)
# and sort the rest according to decreasing normed variance
hvg_rank_sorted_idxs = np.argsort(hvg_rank)
# 04a. Get the sorted rank idxs (the locations of non-zero
# HVG ranks, which are the non-nan ranks, in sorted order)
hvg_rank_sorted_idxs_top_hvgs = [
        _ for _ in
        hvg_rank_sorted_idxs
        if hvg_rank[_]==hvg_rank[_]
        ]
hvg_rank_sorted_gene_names = [
        gene_names[_]
        for _ in
        hvg_rank_sorted_idxs_top_hvgs
        ]
hvg_rank_sorted_var_norm = [
        hvg_var_norm[_]
        for _ in
        hvg_rank_sorted_idxs_top_hvgs
        ]

# 04b. Get the idxs for the genes not in the
# top n HVG values (since these ranks are all
# nan values, they are not necessarily in
# sorted order)
non_hvg_idxs = [
        _ for _ in
        hvg_rank_sorted_idxs
        if hvg_rank[_]!=hvg_rank[_]
        ]
# 04b.i. Get the sorted version of these
# non-HVG indexes by sorting the corresponding
# genes in order of descending normed variance
non_hvg_gene_names = [
        gene_names[_]
        for _ in
        non_hvg_idxs
        ]
non_hvg_var_norm = [
        hvg_var_norm[_]
        for _ in
        non_hvg_idxs
        ]
non_hvg_rank_sorted_idxs_ascending = np.argsort(non_hvg_var_norm)
# 04b.ii. Reverse the order of the idxs
# produced by argsort to put normed variances
# in descending order
non_hvg_rank_sorted_idxs = non_hvg_rank_sorted_idxs_ascending[::-1]
non_hvg_rank_sorted_gene_names = [
        non_hvg_gene_names[_]
        for _ in
        non_hvg_rank_sorted_idxs
        ]
non_hvg_rank_sorted_var_norm = [
        non_hvg_var_norm[_]
        for _ in
        non_hvg_rank_sorted_idxs
        ]

# 04c. Append the non-HVG, rank-ordered
# gene names and normed variances to the HVG,
# rank-ordered gene names and normed
# variances, respectively
ranked_gene_names = hvg_rank_sorted_gene_names + non_hvg_rank_sorted_gene_names
ranked_norm_vars = hvg_rank_sorted_var_norm + non_hvg_rank_sorted_var_norm

# 04d. Print gene names and normed variances in rank order 
ranked_hvg_full_out_fn = f'{hvg_mg_out_dir}/ranked_HVG_list.csv'
ranked_hvg_df = pd.DataFrame(
        data={
            'ranked_HVG_name': ranked_gene_names,
            'ranked_HVG_norm_var': ranked_norm_vars
            }
        )
ranked_hvg_df.to_csv(
        ranked_hvg_full_out_fn,
        index=False
        )

# 05. Pull lists of marker genes
# 05a. Get the list of unique positive marker
# genes in the marker gene database
mg_list_raw = mg_db['marker_genes'].values.tolist()
mg_list = list(
        set(
            [
                item.strip()
                for _ in mg_list_raw
                for item in _.split(',')
                ]
            )
        )
# 05b. Pull negative marker genes in case you want
# to look at these measures with negative marker
# genes included
neg_mg_list_raw = mg_db['negative_markers'].values.tolist()
neg_mg_list_raw = [
        _ for _ in
        neg_mg_list_raw
        if _==_
        ]
neg_mg_list = list(
        set(
            [
                item.strip()
                for _ in neg_mg_list_raw
                for item in _.split(',')
                ]
            )
        )

# 05c. Add the list of unique negative
# marker genes to unique positive marker genes
# in case you decide to use them
#mg_list = list(set(mg_list).union(set(neg_mg_list)))

# 06. For a range of HVG values (in increments
# of 1000), find the number of marker genes in
# the HVG list
n_hvgs_to_check = np.arange(0,37000,1000)
frac_mgs_dict = {}
for n_hvgs_curr in n_hvgs_to_check:
    hvgs_curr = ranked_gene_names[:n_hvgs_curr]
    mgs_in_hvgs_curr = [
            _ for _ in
            mg_list
            if _ in
            hvgs_curr
            ]
    frac_mgs_in_hvgs_curr = 1.0*len(mgs_in_hvgs_curr)/len(mg_list)
    frac_mgs_dict[n_hvgs_curr] = frac_mgs_in_hvgs_curr

# 07. Plot the rank-ordered normed variance
# and mark the positions of 2000 HVGs and 20,000 HVGs,
# superimposed with the fraction of marker genes 
# retained from the list we used
fig,ax = plt.subplots()
# 07a. Plot normed variance and HVG markers
plt.plot(
        np.arange(1,len(ranked_norm_vars)+1),
        ranked_norm_vars,
        color='gray',
        #width=1
        )
pos1 = 2000
pos2 = 20000
ax.plot(
        [pos1]*2,
        [0,np.max(ranked_norm_vars)],
        '-',
        color='dodgerblue',
        linewidth=1,
        label=f'{pos1}'
        )
ax.plot(
        [pos2]*2,
        [0,np.max(ranked_norm_vars)],
        '-',
        color='royalblue',
        linewidth=1,
        label=f'{pos2}'
        )
ax.plot(
        [0,37000],
        [0]*2,
        '--',
        color='lightgray',
        linewidth=1
        )
ax.set_xlabel('HVG Rank')
ax.set_ylabel('Normed Expression Variance')
# 07b. Plot fraction of cell type marker
# genes in HVG list
ax2 = ax.twinx()
ax2.plot(
        [k for k in frac_mgs_dict.keys()],
        [v for v in frac_mgs_dict.values()],
        '--',
        color='darkred',
        linewidth=1,
        )
ax2.tick_params(
        axis='y',
        labelcolor='darkred',
        )
ax2.set_ylabel(
        'Frac. Marker Genes Retained',
        color='darkred',
        rotation=270,
        labelpad=12
        )
ax.legend(loc='center right')
plt_output_full_fn = f'{hvg_mg_out_dir}/ranked_hvg_vars_w_n_hvg_markers.png'
plt.savefig(
        plt_output_full_fn,
        dpi=300
        )
plt.close()

# 08. Plot the expression across
# cell type clusters of certain
# inflammatory and cell type marker
# genes along with their ranked
# positions in the HVG list
infl_gene_list = [
        'MEGF11',
        'CD44',
        'CRP',
        'PPP2R2B',
        'PARP1',
        'FUT4',
        'FUT8',
        'IL1B',
        'IL1RL1',
        'IL7R',
        'IL1R2',
        'IL6',
        'IL1A',
        'IL2RA',
        'TLR2',
        'TLR1',
        'TLR4',
        'TLR5',
        'SEMA6B',
        'SEMA3B',
        'SEMA4B',
        'SEMA4A',
        'SEMA4D',
        'SEMA6A'
        ]
ct_mg_list = [
        'OLIG1',
        'OLIG2',
        'SLC6A3',
        'TH',
        'GFAP',
        ]

# 08a. Define a helper function to
# determine the position of a specified
# gene in an HVG list
def get_hvg_rank(
        gene,
        hvg_ranked_list
        ):
    gene_rank_list = [
            i for i,_ in
            enumerate(hvg_ranked_list)
            if _==gene
            ]
    print(f'preview of hvg ranked list: {hvg_ranked_list[:10]}')
    print(f'current gene rank list: {gene_rank_list}')
    print('\n\n')
    gene_rank = None
    if len(gene_rank_list) > 0:
        gene_rank = gene_rank_list[0]
    return gene_rank

# 08b. Define a helper function to check
# whether each gene in a user-specified list
# is in the HVG list, then if so, to plot
# a 2D UMAP showing the expression for that
# gene across nuclei clusters
def plot_gene_w_rank(
        gene_list,
        hvg_ranked_list,
        adata
        ):
    for gene in gene_list:
        gene_rank = get_hvg_rank(
                gene,
                hvg_ranked_list
                )
        print(f'Current gene: {gene}\nCurrent gene rank: {gene_rank}')
        if gene_rank is not None:
            if gene_rank < N_HVGS:
                fig,ax = plt.subplots()
                sc.pl.embedding(
                        adata,
                        basis='umap',
                        color=gene_prefix+gene,
                        vmin='p0',
                        vmax='p95',
                        color_map='Oranges',
                        size=7.5,
                        legend_loc=None,
                        title=f'{gene}',
                        show=False,
                        ax=ax
                        )
                plt_full_fn = f'{hvg_mg_out_dir}/2D_UMAP_{gene}_HVG_rank_{gene_rank}.png'
                plt.savefig(
                        plt_full_fn,
                        dpi=300
                        )
                plt.close('all')


print(f'{len(ranked_gene_names)}')
# 08c. Inspect HVG ranks for inflammation-related
# genes
plot_gene_w_rank(
        infl_gene_list,
        ranked_gene_names,
        ahc_ct
        )

# 08d. Inspect HVG ranks for cell type marker
# genes
plot_gene_w_rank(
        ct_mg_list,
        ranked_gene_names,
        ahc_ct
        )

                

sys.exit()




























