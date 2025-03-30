import numpy as np
import pandas as pd
import anndata as ad
import scrublet as scr
import scanpy as sc
import scipy.sparse as sp
import matplotlib.pyplot as plt
import os
import sys
import argparse
import yaml
from utils.consts import *

# This script takes CellBender-cleaned, demuxlet-filtered
# expression data, and performs on each individual pooled-sample
# dataset an additional round of doublet detection using Scrublet.
# In this script, each barcode per pooled dataset is assigned a doublet
# score, and the scored expression data is written to an output file.

# Scrublet documentation and code can be found at
# https://www.sciencedirect.com/science/article/pii/S2405471218304745
# https://github.com/swolock/scrublet

# Script setup

# 00. Create argparse object for reading input arguments
parser = argparse.ArgumentParser(description='Doublet score assignment for expression data.')
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
# 01a.iii. Raw input expression data file name
input_data_fn = cfg.get('aggregated_raw_expr_data_fn') 
# 01a.iv. Doublet-scored data output file name
qc_output_fn = cfg.get('dbl_scored_agg_raw_expr_data_fn')
# 01b. Random number generator information
# Note: We are considering this doublet removal step to be
# a separate process, and to keep QC and analysis results
# the same as not using Scrublet when Scrublet detects
# no doublets for removal, we use a separate random state
# here (initialized with a different random seed).
random_seed_init = cfg.get('rand_seed_dbl_rm') 
# 01c. Scrublet information
# 01c.i. Doublet simulation parameters
simulated_doublet_ratio = cfg.get('doublet_simulation_rate')
# 01c.ii. Doublet results output file directory
dbl_results_dir = preprocessed_data_dir + f'/' + cfg.get('scrublet_results_dir')
# 01c.ii.A. Append the random seed to the results dir
dbl_results_dir += f'{random_seed_init}'
if not os.path.isdir(dbl_results_dir):
    os.system(f'mkdir -p {dbl_results_dir}')
# 01c.iii. Doublet score threshold file name
dbl_threshold_out_fn = cfg.get('doublet_threshold_fn')

# 02. Import expression data
ahc_fn = f'{preprocessed_data_dir}/{input_data_fn}'
adatas_human = ad.read_h5ad(ahc_fn)

# 02a. Make a working copy of the barcode
# annotation DataFrame for doublet score assignment
ahc_obs_cp = adatas_human.obs.copy()

# 03. Import per-pool demuxlet multi-individual
# doublet counts to use for single-individual
# doublet count estimation
mi_doublet_fn = f'{preprocessed_data_dir}/demuxlet_multi_ind_doublet_counts.csv'
mi_doublet_df = pd.read_csv(mi_doublet_fn,
        index_col=0)
mi_doublet_df['pool_name'] = mi_doublet_df['Pool_ID'].values.tolist()

# 04. Initialize a DataFrame to hold the Scrublet-generated
# doublet detection thresholds for each patient in each pool.
# This information will be adjusted after manual inspection
# of doublet scores for simulated doublets and observed
# transcriptomes, and the adjusted values will be used
# to call any possible remaining single-individual doublets
# for each patient in a downstream processing step.
doublet_threshold_df = pd.DataFrame(
        columns=[
            'patient_ID',
            'pool_name',
            'scrublet_generated_doublet_threshold',
            'manually_adjusted_doublet_threshold'
            ]
        )
dt_idx = 0

# 05. Initialize the random state with the input
# random seed. This generator will be used to generate
# random seeds for each call to Scrublet.
numpy_random_state = np.random.RandomState(
        seed=random_seed_init
        )
# 05a. Define a helper function to get a random
# seed from the current random number generator
def get_rand_seed(random_number_generator):
    rand_seed_curr = numpy_random_state.randint(
            0,
            INT32_MAX
            )
    return rand_seed_curr

# 06. Perform doublet detection on each pool in the expression
# data
# 06a. Split expression data by pool and patient so that
# each individual sample is processed on its own, per
# recommendations in Scrublet documentation (because
# doublets are detected through transcriptomic simulations
# based on actual data, and unaccounted-for interindividual
# or inter-pool differences may make these simulations
# less accurate).
# 06a.i. Pull list of pooled dataset names
pools = list(
        np.unique(
            adatas_human.obs['pool_name'].values.tolist()
            )
        )
# 06a.ii. Process data for each individual in each pool
for pool in pools:
    print(f'Running scrublet on pool {pool}...')
    ahc_curr = adatas_human[
            adatas_human.obs['pool_name'] == pool,
            :
            ].copy()
    # 06b. Calculate the expected number of single-individual
    # doublets for the entire pool, then compute the expected
    # number per individual under the assumption that
    # that single-individual (si) doublets are uniformly likely
    # across individuals (since this is a technical artifact)
    n_mi_doublets = mi_doublet_df.loc[
            mi_doublet_df['pool_name'] == pool
            ]['N_DBL_expected_patients'].values.tolist()[0]
    print(f'\tNumber of multi-individual doublets from demuxlet: {n_mi_doublets}')
    pts_pool_curr = list(
            np.unique(
                ahc_curr.obs['patient_ID'].values.tolist()
                )
            )
    n_pts_pool_curr = len(pts_pool_curr)
    print(f'\tNumber of patients in this pool: {n_pts_pool_curr}')
    
    # 06b.i. Use the multi-individual doublet number to estimate
    # the total expected number of si doublets in the pool.
    # Here, we assume that multi-individual doublets from any 
    # combination of individuals is equally likely in our derivation.  
    # Note: for pools containing tissue from only one individual, the number
    # of multi-individual doublets has to be zero by definition. In these
    # cases, we are unable to estimate the number of si doublets and will
    # use Scrublet's default parameter.
    n_si_doublets = None
    if n_pts_pool_curr > 1:
        n_si_doublets = (1/(n_pts_pool_curr-1))*n_mi_doublets
        print(f'\tEstimating the number of single-individual doublets in this multi-individual pool.' + \
                f'\n\t\tEstimated number of single-individual doublets per pool: {n_si_doublets:0.2f}')   
        # 06b.ii. Estimate the expected number of si doublets per individual
        n_si_doublets_per_patient = n_si_doublets/n_pts_pool_curr
        print(f'\t\tEstimated number of single-individual doublets per patient: {n_si_doublets_per_patient:0.2f}')

    # 06c. Estimate doublets for individuals in the current pool
    for pt_curr in pts_pool_curr:
        print(f'\n\t\tProcessing patient {pt_curr} in this pool...')
        # 06c.i. Get expression data for the current donor
        ahc_pt_curr = ahc_curr[
                ahc_curr.obs['patient_ID'] == pt_curr,
                :
                ]
        ahc_pt_curr_x = ahc_pt_curr.X.copy()

        # 06c.ii. Compute the expected si doublet fraction if possible
        # and use it to initialize the Scrublet object.

        # From our intial inspections, we can see that the default simulated
        # doublet rate (2) is not generating doublets with very high scores, even
        # though most samples have most cell types represented. This may suggest
        # that the default rate is not quite high enough to cover all possible
        # doublet configurations (possibly just based on the imbalance of cell type
        # proportions per patient).
        # So, here, we will initialize sim_doublet_ratio to a slightly higher
        # value to try and make sure the least-represented cell types are also 
        # captured.
        if n_si_doublets:
            n_barcodes_tot_pt_curr = ahc_pt_curr_x.shape[0]
            frac_si_doublets_pt_curr = n_si_doublets_per_patient*1.0/n_barcodes_tot_pt_curr
            print(f'\t\tEstimated fraction of single-individual doublets for this patient: {frac_si_doublets_pt_curr:0.2f}')
            print(f'\t\tInitializing Scrublet object...')
            scrub = scr.Scrublet(
                    ahc_pt_curr_x,
                    expected_doublet_rate=frac_si_doublets_pt_curr,
                    sim_doublet_ratio=simulated_doublet_ratio,
                    random_state=get_rand_seed(numpy_random_state)
                    )
        else:
            print(f'\t\tInitializing Scrublet object...')
            scrub = scr.Scrublet(
                    ahc_pt_curr_x,
                    sim_doublet_ratio=simulated_doublet_ratio,
                    random_state=get_rand_seed(numpy_random_state) # rand_seed_curr
                    )
            
        # 06d. Run the default Scrublet doublet detection pipeline
        # Note: because we had to keep a larger number of HVGs to 
        # see certain marker genes in this dataset, we will also retain
        # a larger-than-default number of genes for doublet detection, by
        # using a lower-than-default minimum threshold on gene variability
        # (the default is 85)
        print(f'\t\tPerforming Scrublet workflow...')
        doublet_scores,predicted_doublets = scrub.scrub_doublets(
                min_counts=2,
                min_cells=3,
                min_gene_variability_pctl=75,
                n_prin_comps=30
                )

        # 06e. Store the Scrublet-generated doublet detection threshold
        # for the current donor
        doublet_threshold_df.at[dt_idx,'patient_ID'] = pt_curr
        doublet_threshold_df.at[dt_idx,'pool_name'] = pool
        doublet_threshold_df.at[dt_idx,'scrublet_generated_doublet_threshold'] = scrub.threshold_.astype('float32')
        dt_idx += 1

        # 06f. Add doublet information to the associated barcodes
        for o_idx_i,o_idx in enumerate(ahc_pt_curr.obs.index.values.tolist()):
            ahc_obs_cp.at[o_idx,'doublet_score'] = doublet_scores[o_idx_i]
            ahc_obs_cp.at[o_idx,'is_doublet'] = (
                    0 if predicted_doublets[o_idx_i]=='False'
                    else 1)

        # 06g. Generate visualizations as performance checks
        # to help adjust the output detection threshold
        # for each patient in downstream filtering (to make
        # sure the thresholds for doublet calling appear to
        # actually remove doublets)

        # 06g.i. Make an output directory for Scrublet's plotting
        # functions
        scr_vis_output_dir = f'{dbl_results_dir}/scrublet_plots'
        if not os.path.isdir(scr_vis_output_dir):
            os.system(f'mkdir -p {scr_vis_output_dir}')
        scr_vis_hist_output_dir = f'{scr_vis_output_dir}/doublet_score_histograms'
        if not os.path.isdir(scr_vis_hist_output_dir):
            os.system(f'mkdir -p {scr_vis_hist_output_dir}')
        # 06g.ii. Plot histograms of the simulated doublet
        # and transcriptome doublet scores
        fig,axs = scrub.plot_histogram()
        hist_fn = f'{scr_vis_hist_output_dir}/doublet_scores__pool_{pool}_pt_{pt_curr}.png'
        plt.savefig(
                hist_fn,
                dpi=300
                )
        del fig,axs
        # 06g.iii. Plot the 2D UMAP for the current donor and pool,
        # with labels showing predicted doublets
        scr_vis_umap_output_dir = f'{scr_vis_output_dir}/2D_UMAPs'
        if not os.path.isdir(scr_vis_umap_output_dir):
            os.system(f'mkdir -p {scr_vis_umap_output_dir}')
        scrub.set_embedding(
                'UMAP',
                scr.get_umap(
                    scrub.manifold_obs_,
                    10,
                    min_dist=0.3
                    )
                )
        fig,axs = scrub.plot_embedding(
                'UMAP',
                order_points=True
                )
        umap_fn = f'{scr_vis_umap_output_dir}/doublet_score_umaps__pool_{pool}_pt_{pt_curr}.png'
        plt.savefig(
                umap_fn,
                dpi=300
                )
        del fig,axs
    
    print('\n')

# 07. Add observation annotations with doublet scores back into expression data
adatas_human.obs = ahc_obs_cp.copy()
del ahc_obs_cp

# 08. Write scored expression data and doublet threshold DataFrame to the
# output directory
# 08a. Write expression data to file
ahc_out_fn = f'{preprocessed_data_dir}/{qc_output_fn}'
adatas_human.write(ahc_out_fn)

# 08b. Write doublet threshold DataFrame to file
doublet_threshold_fn = f'{dbl_results_dir}/{dbl_threshold_out_fn}'
doublet_threshold_df.to_csv(doublet_threshold_fn,index=True)

sys.exit()

