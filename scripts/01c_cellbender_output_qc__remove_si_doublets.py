import numpy as np
import pandas as pd
import anndata as ad
import scrublet as scr
import matplotlib.pyplot as plt
import os
import sys
import argparse
import yaml
from utils.consts import *
from utils import reformat_pt_ids as rpi

# This script takes adjusted doublet thresholds for each
# sample (i.e. each patient in each unique pool) for specified
# expression data and uses it to detect and remove barcodes that
# resemble doublets (according to the results of Scrublet).
# It saves this filtered expression data to be used for downstream
# processing.


# Script setup

# 00. Create argparse object for reading input arguments
parser = argparse.ArgumentParser(description='Double removal for aggregated, raw expression data.')
parser.add_argument(f'--config-yaml-path',type=str,required=True)
args = parser.parse_args()

# 01. Get analysis parameters from configuration file
cfg = None
with open(args.config_yaml_path) as f:
    cfg = yaml.safe_load(f)
# 01a. Dataset information
# 01a.i. Data root directory
qc_root_dir = cfg.get('root_dir')
# 01a.ii. Data directory containing partially preprocessed expression data files
preprocessed_data_dir = qc_root_dir + f'/' + cfg.get('preprocessed_data_dir')
# 01a.iv. Doublet-scored input data file name
preprocessed_data_fn = cfg.get('dbl_scored_agg_raw_expr_data_fn')
# 01a.v. Doublet-removed output data file name
qc_output_fn = cfg.get('dbl_rm_agg_raw_expr_data_fn')

# 01b. Doublet scoring information
# 01b.i. Directory with Scrublet results
scrublet_results_dir = preprocessed_data_dir + f'/' + cfg.get('scrublet_results_dir')
# 01b.i.A. Append the Scrublet random seed to directory name
rand_seed_dbl = cfg.get('rand_seed_dbl_rm')
scrublet_results_dir += f'{rand_seed_dbl}'
# 01b.ii. Directory with manually validated Scrublet threshold file
validated_doublet_threshold_dir = qc_root_dir + f'/' + cfg.get('validated_doublet_threshold_dir')
# 01b.iii. Validated doublet threshold file fn
validated_doublet_threshold_fn = cfg.get('validated_doublet_threshold_fn')
# 01b.iv. Name of output file with information about removed doublets
removed_dbl_output_fn = cfg.get('removed_doublet_fn')

# Doublet removal

# 02. Import data
# 02a. Import expression data
ahc_fn = f'{preprocessed_data_dir}/{preprocessed_data_fn}'
adatas_human = ad.read_h5ad(ahc_fn)

# 02b. Make a working copy of the
# barcode annotation DataFrame
ahc_obs_cp = adatas_human.obs.copy()
# 02b.i. Drop the 'is_doublet' column,
# which flags barcodes based on non-validated,
# Scrublet-selected thresholds, to avoid confusion
ahc_obs_cp.drop(columns=['is_doublet'],
        inplace=True)
# 02b.ii. Reformat the patient ID list as well
ahc_obs_cp['patient_ID'] = ahc_obs_cp['patient_ID'].astype('object')
ahc_obs_cp['patient_ID'] = [f'{_}' for _ in ahc_obs_cp['patient_ID'].values.tolist()]
ahc_obs_cp['patient_ID'] = rpi.reformat_pt_ids(ahc_obs_cp['patient_ID'])

# 02c. Import validated doublet threshold data
dt_full_fn = f'{validated_doublet_threshold_dir}/{validated_doublet_threshold_fn}'
doublet_threshold_df = pd.read_csv(
        dt_full_fn,
        index_col=0
        )

# 02d. Initialize a list of barcodes flagged
# as single-individual doublets
doublet_barcodes = []

# 02e. Initialize DataFrame to hold the number
# of doublets removed for each sample
doublets_to_remove_df = pd.DataFrame(
        columns = [
            'patient_ID',
            'pool_name',
            'n_doublets_removed'
            ]
        )
dr_idx = 0

# 02f. Reformat the patient IDs in the doublet threshold
# DataFrame so patients can be matched
doublet_threshold_df['patient_ID'] = doublet_threshold_df['patient_ID'].astype('object')
doublet_threshold_df['patient_ID'] = [
        f'{_}' for _ in doublet_threshold_df['patient_ID'].values.tolist()
        ]
doublet_threshold_df['patient_ID'] = rpi.reformat_pt_ids(
        doublet_threshold_df['patient_ID']
        )

# 03. Apply doublet thresholds to each sample's barcodes
pools = list(
        np.unique(
            ahc_obs_cp['pool_name'].values.tolist()
                )
            )
# 03a. Iterate over pools
for pool in pools:
    pts_curr = list(
            np.unique(
                ahc_obs_cp.loc[
                    ahc_obs_cp['pool_name'] == pool
                    ]['patient_ID'].values.tolist()
                )
            )
    # 03b. Iterate over patients
    for pt in pts_curr:
        # 03c. Pull barcode annotation
        # DataFrame for the current sample
        ahc_curr = ahc_obs_cp.loc[
                (ahc_obs_cp['pool_name'] == pool)
                &
                (ahc_obs_cp['patient_ID'] == pt)
                ]
    
        # 03d. Get the adjusted doublet threshold 
        # for the current sample
        adj_doublet_threshold = doublet_threshold_df.loc[
                (doublet_threshold_df['pool_name'] == pool)
                &
                (doublet_threshold_df['patient_ID'] == pt)
                ]['manually_adjusted_doublet_threshold'].values.tolist()[0]
        
        # 03e. Find the barcodes for which the doublet score
        # is greater than the adjusted threshold
        doublet_barcodes_pt_curr = ahc_curr.loc[
                ahc_curr['doublet_score'] > adj_doublet_threshold
                ].index.values.tolist()

        # 03f. Add these barcodes to the list of doublets
        # to remove
        doublet_barcodes.extend(doublet_barcodes_pt_curr)

        # 03g. Add the number and fraction of removed doublets for this
        # patient in this pool to the doublet removal DataFrame
        doublets_to_remove_df.at[dr_idx,'patient_ID'] = pt
        doublets_to_remove_df.at[dr_idx,'pool_name'] = pool
        doublets_to_remove_df.at[dr_idx,'n_doublets_removed'] = len(doublet_barcodes_pt_curr)
        doublets_to_remove_df.at[dr_idx,'frac_doublets_removed'] = len(doublet_barcodes_pt_curr)*1.0/len(ahc_curr)
        dr_idx += 1

# 04. Remove doublets from expression data and reformat
# it for downstream analysis
# 04a. Get the list of barcodes to keep
barcodes_to_keep = [
        _ for _ in
        ahc_obs_cp.index.values.tolist()
        if _ not in doublet_barcodes
        ]
# 04b. Delete the working copy of the observation annotation DataFrame
# to save space
del ahc_obs_cp

# 04c. Make a copy of the expression data with only barcodes that
# were not called as doublets
adatas_human_dbl_filt = adatas_human[
        barcodes_to_keep,
        :
        ].copy()

# 04d. Remove the columns with doublet score information
# from the filtered expression data
cols_to_keep = [
        _
        for _ in adatas_human_dbl_filt.obs.columns.values.tolist()
        if 'doublet' not in _]
adatas_human_dbl_filt.obs = adatas_human_dbl_filt.obs[cols_to_keep].copy()

# 05. Save the filtered expression data and doublet removal information
# 05a. Save expression data
ahc_out_fn = f'{preprocessed_data_dir}/{qc_output_fn}'
adatas_human_dbl_filt.write(ahc_out_fn)

# 05b. Write removed doublet DataFrame to file
doublets_to_remove_fn = f'{scrublet_results_dir}/{removed_dbl_output_fn}'
doublets_to_remove_df.to_csv(doublets_to_remove_fn,index=True)

# 06. Generate plots to inspect filtering results
# 06a. Plot a histogram of the number of barcodes removed per sample
plt.figure()
binsize = 50
n_dbl_rem = doublets_to_remove_df['n_doublets_removed'].values.tolist()
max_n_rem = np.max(n_dbl_rem)
plt.hist(
        n_dbl_rem,
        bins=np.arange(0,max_n_rem+2*binsize,binsize),
        alpha=0.8
        )
plt.xlabel('N doublets removed')
plt.ylabel('Frequency')
n_dbl_rem_fn = f'{scrublet_results_dir}/hist__n_doublets_removed_per_sample.png'
plt.savefig(
        n_dbl_rem_fn,
        dpi=300
        )
plt.close()
del n_dbl_rem,max_n_rem,binsize

# 06b. Plot a histogram of the fraction of barcodes removed per sample
plt.figure()
frac_dbl_rem = doublets_to_remove_df['frac_doublets_removed'].values.tolist()
binsize = 0.01
plt.hist(
        frac_dbl_rem,
        bins=np.arange(0,1.0+2*binsize,binsize),
        alpha=0.8
        )
plt.xlabel('Frac. doublets removed')
plt.ylabel('Frequency')
frac_dbl_rem_fn = f'{scrublet_results_dir}/hist__frac_doublets_removed_per_sample.png'
plt.savefig(
        frac_dbl_rem_fn,
        dpi=300
        )
plt.close()



