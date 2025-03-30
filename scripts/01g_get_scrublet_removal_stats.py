import numpy as np
import pandas as pd
import os
import sys
import argparse
import yaml

# This script computes summary statistics
# about the number of single-individual
# doublets removed from expression data
# with Scrublet.

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
# 01a.ii. Dataset parent directory
dataset_parent_dir = cfg.get('data_parent_dir')
# 01a.iii. Data directory containing preprocessed expression data files
preprocessed_data_dir = qc_root_dir + f'/' + dataset_parent_dir + f'/' + cfg.get('preprocessed_data_dir')
# 01a.iv. Doublet-scored input data file name
preprocessed_data_fn = cfg.get('dbl_scored_agg_raw_expr_data_fn')
# 01a.v. Doublet-removed output data file name
qc_output_fn = cfg.get('dbl_rm_agg_raw_expr_data_fn')

# 01b. Doublet scoring information
# 01b.i. Directory for Scrublet-generated
# results
scrublet_results_dir = preprocessed_data_dir + f'/' + cfg.get('scrublet_results_dir')
# 01b.i.A. Append the Scrublet random seed
# to complete the directory name
scrublet_random_seed = cfg.get('rand_seed_dbl_rm')
scrublet_results_dir += f'{scrublet_random_seed}'
# 01b.ii. Directory containing validated
# Scrublet detection threshold file
# (which may differ from the output
# directory above)
validated_doublet_threshold_dir = qc_root_dir + f'/' + dataset_parent_dir + f'/' + cfg.get('validated_doublet_threshold_dir')
# 01b.iii. Validated doublet threshold file fn
validated_doublet_threshold_fn = cfg.get('validated_doublet_threshold_fn')
# 01b.iv. Name of output file with information about removed doublets
removed_dbl_output_fn = cfg.get('removed_doublet_fn')
# 01b.v. Directory containing publication run
# information about removed doublets
pub_run_removed_dbl_output_dir = qc_root_dir + f'/' + dataset_parent_dir + f'/' + cfg.get('publication_run_removed_doublet_dir') 
# 01b.vi. Publication run file with information
# about removed doublets
pub_run_removed_dbl_output_fn = cfg.get('publication_run_removed_doublet_fn')

# 01c. Specify a helper function for computing
# stats
def get_stats(vals):
    """
    This function computes
    a set of statistics on
    a user-specified list of values
    and returns them in a dictionary

    vals: a list of numeric values
    stat_dict: the dictionary containing
    summary statistics for the input value list
    """

    stat_dict = {}
    stat_dict['min'] = np.min(vals)
    stat_dict['max'] = np.max(vals)
    stat_dict['mean'] = np.mean(vals)
    stat_dict['median'] = np.median(vals)
    pctile_25 = np.percentile(vals,25)
    pctile_75 = np.percentile(vals,75)
    stat_dict['IQR'] = f'({pctile_25},{pctile_75})'

    return stat_dict


# 02. Import Scrublet removal information
# from the current run and the publication
# run, and for each, produce summary information
removal_dbl_full_fn_dict = {
        'curr_run': {
            fn_tag: '',
            full_fn: f'{scrublet_results_dir}/{removed_dbl_output_fn}'
            },
        'pub_run': {
            fn_tag: '__AW_pub_run',
            full_fn; f'{pub_run_removed_dbl_output_dir}/{pub_run_removed_dbl_output_fn}'
            }
        }
for file_type, file_info_dict in removal_dbl_full_fn_dict.items():
    fn_tag = file_info_dict['fn_tag']
    removal_dbl_full_fn = file_info_dict['full_fn']
    removal_df = None
    ndb_stats_df = pd.DataFrame(
            columns=['N_Dbl_Removed','Frac_Dbl_Removed']
            )
    if os.path.isfile(removal_dbl_full_fn):
        removal_df = pd.read_csv(
                removal_dbl_full_fn,
                index_col=0
                )
        # 02a. Compute the minimum, maximum,
        # mean, median, and IQR for the number
        # of doublets removed by Scrublet
        # for the nuclei assigned to each
        # individual for each library
        n_dbl = removal_df['n_doublets_removed'].values.tolist()
        frac_dbl = removal_df['frac_doublets_removed'].values.tolist()
        
        n_stats_dict = get_stats(n_dbl)
        frac_stats_dict = get_stats(frac_dbl)
        
        # 02b. Store the summary statistics in a DataFrame
        ndb_stats_df.at['min','N_Dbl_Removed'] = n_stats_dict['min']
        ndb_stats_df.at['max','N_Dbl_Removed'] = n_stats_dict['max']
        ndb_stats_df.at['mean','N_Dbl_Removed'] = n_stats_dict['mean']
        ndb_stats_df.at['median','N_Dbl_Removed'] = n_stats_dict['median']
        ndb_stats_df.at['IQR','N_Dbl_Removed'] = n_stats_dict['IQR']
        
        ndb_stats_df.at['min','Frac_Dbl_Removed'] = frac_stats_dict['min']
        ndb_stats_df.at['max','Frac_Dbl_Removed'] = frac_stats_dict['max']
        ndb_stats_df.at['mean','Frac_Dbl_Removed'] = frac_stats_dict['mean']
        ndb_stats_df.at['median','Frac_Dbl_Removed'] = frac_stats_dict['median']
        ndb_stats_df.at['IQR','Frac_Dbl_Removed'] = frac_stats_dict['IQR']
        
        # 02c. Save the summary statistic DataFrame to file
        ndb_stats_full_fn = f'{preprocessed_data_dir}/{scrublet_results_dir}/rem_dbl_summ_stats_ppt_ppool{fn_tag}.csv'
        ndb_stats_df.to_csv(
                ndb_stats_full_fn
                )

sys.exit()

