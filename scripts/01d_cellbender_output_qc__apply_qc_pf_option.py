import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import scipy.sparse as sp
import matplotlib.pyplot as plt
import os
import harmonypy as hm
import sys
import pickle
import argparse
import yaml
from utils.consts import *
from utils.get_rand_seed import get_rand_seed
from utils.sparse_funcs import sparse_mean_var_major_axis, sparse_mean_var_minor_axis, sparse_mean_variance_axis, get_mean_var__dtype_mod, scale_array__dtype_mod 
from utils.assess_qc_metrics import compute_qc_metrics, plot_qc_metric_distributions

# This script takes snRNA-seq data that has been aggregated
# and undergone cleaning and doublet-removal steps
# (CellBender, demuxlet, and Scrublet) and applies
# additional quality control (QC) and preprocessing steps
# to it, including
# (1) Filtering to remove nuclei of potentially poor quality 
#     based on certain heuristics (fraction of mitochondrial
#     reads, number of genes, and number of unique molecular
#     identifiers [UMIs] per nucleus)
# (2) Dimensionality reduction by highly variable gene (HVG)
#     selection and principal component analysis (PCA)
# (3) Batch correction across sequencing libraries
#     using Harmony
# (4) Propagation of dimensionality reduced expression data
#     for downstream analysis
# (5) Nearest-neighbor computation, UMAP dimensionality reduction,
#     and initial Leiden cluster computation
# This script then saves the preprocessed expression data
# for final downstream adjustments and rescaling.

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
output_data_batch_corr_fn = cfg.get('preproc_exp_data_batch_corr_scaled_fn')
output_data_pf_only_fn = cfg.get('preproc_exp_data_filt_pf_only_fn')

# 01b. Preprocessing parameters
# 01b.i. Patient IDs to remove
# (e.g. because of below-power-threshold nuclear yield)
pt_IDs_to_remove = cfg.get('patients_to_remove')
# 01b.iii. Harmony batch correction variable
harmony_batch_var = cfg.get('harmony_batch_corr_var')
# 01b.iv. Whether to filter based on mitochondrial
# read fraction
FILTER_ON_MITO_READS = cfg.get('mito_filtering')
# 01b.v. Whether to only perform QC heuristic filtering
# 'pass/fail' assignment to every barcode, save this
# information, and then exit instead of performing
# the full QC run
PASS_FAIL_QC_ONLY = cfg.get('pass_fail_only')
# 01b.vi. Dictionary with QC hyperparameters
qc_hyperparam_dict = cfg.get('qc_hyperparams')
# 01b.vii. Number of HVGs to retain
N_HVGS = qc_hyperparam_dict['n_hvgs']
# 01b.viii. Scale cap on HVG expression values
MAX_VAL_HVG_SCALING = qc_hyperparam_dict['hvg_expression_scale_cap_val']
if MAX_VAL_HVG_SCALING == 'None':
    MAX_VAL_HVG_SCALING = None
# 01b.ix. Number of nearest neighbors to use
# for network computation and cluster generation
N_NEIGHBORS = int(qc_hyperparam_dict['n_neighbors'])

# 01c. Random state setup (random state was seeded
# here and preserved throughout subsequent processing
# scripts where possible)
# 01c.i. Random seed
random_seed_init = cfg.get('rand_seed_qc_preproc')
# 01c.ii. Output directory for RNG state files
random_state_dir = qc_root_dir + f'/' + cfg.get('rng_state_dir')
# 01c.iii. Output RNG state file name
random_state_fn = cfg.get('rng_state_fn__01_preproc')
# 01c.iv. Initialize random state
numpy_random_state = np.random.RandomState(
        seed=random_seed_init
        )
# 01c.v. Set up random state output directory
# if necessary
if not os.path.isdir(random_state_dir):
    os.system(f'mkdir -p {random_state_dir}')

# 01d. Set figure params
sc.set_figure_params(fontsize=10, dpi_save=300, vector_friendly=True)

# Perform QC

# 02. Read expression data
ahc_fn = f'{preprocessed_data_dir}/{input_data_fn}'
adatas_human = ad.read_h5ad(ahc_fn)

# 03. Set the datatype of adatas_human to one with
# <64 bit depth if possible, to free up memory
max_count_val_for_dtype = np.max(adatas_human.X.copy())
working_precision =''
if max_count_val_for_dtype < INT32_MAX:
    working_precision = '32'
else:
    print(f'This count matrix requires 64-bit precision, and many preprocessing steps',
            f' (such as HVG detection) may not function as a result.')
    working_precision = '64'
working_dtype_int = f'uint{working_precision}'
working_dtype_float = f'float{working_precision}'
print(f'Integer dtype for processing: {working_dtype_int}\nFloat dtype for processing: {working_dtype_float}')
adatas_human.X = adatas_human.X.copy().astype(working_dtype_int)

# 04. Perform a standard single-nucleus QC workflow on the aggregated data

# 04a. Define a helper function to compute metrics for QC, including n total molecules, 
# n unique genes, and fraction of mitochondrial reads
# per barcode. This function adds these measures to the input anndata object inplace.
def compute_qc_metrics(adata: ad.AnnData):
    # 04a.i. Calculate a few metrics and look at them
    if 'mito_frac' not in adata.obs.keys():
        sys.stderr.write('Calculating mito_frac, assuming mito genes start with\ngenome prefix plus "mt-" or "MT-"')
        # 04a.ii. This works for our "<genome>________________MT-" format for mitochondrial genes in the
        # human case. Check that it also works for mouse genome data.
        if adata.var_names.str.contains('mt-',case=False).sum() == 0:
            sys.stderr.write('WARNING: no genes found that start with "mt-"')
        adata.obs['mito_frac'] = np.array(
            adata.X[:, adata.var_names.str.contains('mt-',case=False)].sum(axis=1)
        ).squeeze() / np.array(adata.X.sum(axis=1) + 1e-10).squeeze()
    adata.obs['n_umi'] = np.array(adata.X.sum(axis=1)).squeeze()
    adata.obs['n_gene'] = np.array((adata.X > 0).sum(axis=1)).squeeze()

# 04b. Function to perform the full QC workflow (including filtration based on the above QC metrics, 
# followed by population standardization and detection of highly variable genes, then for the subset
# of highly variable genes, standardization of expression per gene across the sample population,
# followed by linear and nonlinear dimensionality reduction and clustering of barcodes based on
# expression profile similarity.
def qc_and_typical_workflow(adata: ad.AnnData, 
                            umi_percentile: float = 95, 
                            gene_percentile: float = 95, 
                            n_gene_min: float = 0,
                            mito_frac_percentile: float = 90,
                            hvg_flavor: str = 'seurat_v3',
                            hvg_scaling_cap: float = None,
                            n_hvgs_to_keep: int = 2000,
                            n_pcs: int = 25,
                            n_neighbors = 20,
                            leiden_resolution: float = 0.5,
                            working_dtype_int: str = 'int64',
                            working_dtype_float: str = 'float64',
                            scaling_accumulator_dtype: str = 'float64',
                            pass_fail_only_mode: bool = False):
    """Do most basic cell QC and a typical scanpy workflow"""
    # 04b.i. Indicate that instances of numpy_random_state
    # refer to the global variable (global-ness is usually
    # implied by variables not defined in a function, unless
    # they are changed in the function, which is why we need
    # to make this declaration here)
    global numpy_random_state
    # 04b.ii. Calculate metrics for QC filtering
    compute_qc_metrics(adata)
    # 04b.iii. Plot distributions of QC metrics
    plot_qc_metric_distributions(adata,
            output_dir = preprocessed_data_dir)
    # 04b.iv. If pass-fail-only mode is set, translate QC metrics
    # to fail flags and return the resulting AnnData object
    if pass_fail_only_mode:
        sys.stderr.write(f'Running QC filtration pass-fail labeling only for data')
        print(f'Data prior to QC: {adata.shape[0]} cells, {adata.shape[1]} genes.')
        adata_obs_cp = adata.obs.copy()
        # 04b.v. Print percentile-based thresholds
        # for each QC metric
        print(f'\tThreshold values for QC metrics:')
        pctile_n_umi = np.percentile(
                adata_obs_cp['n_umi'],
                q=umi_percentile
                )
        print(f'\tN UMI: < {pctile_n_umi}')
        pctile_n_genes = np.percentile(
                adata_obs_cp['n_gene'],
                q=gene_percentile
                )
        print(f'\tN genes: > {n_gene_min}, < {pctile_n_genes}')
        pctile_mito_frac = np.percentile(
                adata_obs_cp['mito_frac'],
                q=mito_frac_percentile
                )
        print(f'\tMito. frac.: < {pctile_mito_frac}\n\n')
        del pctile_n_umi,pctile_n_genes,pctile_mito_frac
        # 04b.v. Perform n_umi pass-fail
        obs_n_umi_pass_idxs = adata_obs_cp[
                adata_obs_cp['n_umi']
                < np.percentile(
                    adata_obs_cp['n_umi'],
                    q=umi_percentile)
                ].index
        print(f'Example index, n_umi: {obs_n_umi_pass_idxs[0]}')
        adata_obs_cp['qc_fail_n_umi_upper_pct_bound'] = 1
        adata_obs_cp.loc[obs_n_umi_pass_idxs,
                'qc_fail_n_umi_upper_pct_bound'] = 0
        # 04b.vi. Perform n_gene pass-fail
        # Flag barcodes with too many total genes
        obs_n_gene_ub_pass_idxs = adata_obs_cp[
                adata_obs_cp['n_gene']
                < np.percentile(
                    adata_obs_cp['n_gene'],
                    q=gene_percentile)
                ].index
        print(f'Example index, n_gene upper bound: {obs_n_gene_ub_pass_idxs[0]}')
        adata_obs_cp['qc_fail_n_gene_upper_pct_bound'] = 1
        adata_obs_cp.loc[obs_n_gene_ub_pass_idxs,
                'qc_fail_n_gene_upper_pct_bound'] = 0
        # 04b.vii. Flag barcodes with too few total genes
        obs_n_gene_lb_pass_idxs = adata_obs_cp[
                adata_obs_cp['n_gene']
                >= n_gene_min
                ].index
        print(f'Example index, n_gene lower bound: {obs_n_gene_lb_pass_idxs[0]}')
        adata_obs_cp['qc_fail_n_gene_lower_n_bound'] = 1
        adata_obs_cp.loc[obs_n_gene_lb_pass_idxs,
                'qc_fail_n_gene_lower_n_bound'] = 0
        # 04b.viii. Perform mito fraction pass-fail
        obs_mito_frac_pass_idxs = adata_obs_cp[
                adata_obs_cp['mito_frac']
                < np.percentile(
                    adata_obs_cp['mito_frac'],
                    q=mito_frac_percentile)
                ].index
        print(f'Example index, mito_frac upper bound: {obs_mito_frac_pass_idxs[0]}')
        adata_obs_cp['qc_fail_mito_frac_upper_pct_bound'] = 1
        adata_obs_cp.loc[obs_mito_frac_pass_idxs,
                'qc_fail_mito_frac_upper_pct_bound'] = 0
        # 04b.ix. Re-assign observation metadata to the AnnData object
        # and return it
        adata_qc = adata.copy()
        adata_qc.obs = adata_obs_cp.copy()
        return adata_qc

    else:
        # 04b.x. Get rid of outlier cells
        sys.stderr.write(f'{adata.shape[0]} cells before cell QC\n')
        adata_qc = adata[(adata.obs['n_umi'] < np.percentile(adata.obs['n_umi'], q=umi_percentile))
                         & (adata.obs['n_gene'] < np.percentile(adata.obs['n_gene'], q=gene_percentile))
                         & (adata.obs['mito_frac'] < 
                            np.percentile(adata.obs['mito_frac'], q=mito_frac_percentile))
                         & (adata.obs['n_gene'] >= n_gene_min)].copy()
        sys.stderr.write(f'{adata_qc.shape[0]} cells after cell QC\n')
        # 04b.xi. Typical scanpy workflow
        # Dataset-wide standardization for highly variable gene (hvg) detection
        # We default to using 'seurat_v3'-flavored hvg detection, since it
        # incorporates some harmonization into its approach.
        # This flavor expects raw counts data (not lognormalized data as the other
        # flavors do).
        # Note: there is a batch key that can be used to perform some batch
        # correction during hvg selection. However, pre-Harmony batch
        # correction measures are not recommended, since they can introduce
        # issues that reduce the overall performance of batch correction, so
        # we do not use these batch keys here.
        adata_qc.X = adata_qc.X.copy()
        sys.stderr.write(f'Detecting highly variable genes with {hvg_flavor}...\n')
        if hvg_flavor != 'seurat_v3':
            sys.stderr.write(f'Lognormalizing aggregated count data\n')
            sc.pp.normalize_total(adata_qc)
            sc.pp.log1p(adata_qc)
        # 04b.xii. Detect hvgs
        print(f'Finding top {n_hvgs_to_keep} highly variable genes...\n')
        sc.pp.highly_variable_genes(adata_qc, 
                                    n_top_genes=n_hvgs_to_keep,
                                    flavor=hvg_flavor)
        adata_qc.obsm['hvgX'] = adata_qc.X[:, adata_qc.var['highly_variable']].copy()
        # 04b.xiii. Do population-wide expression level scaling per hvg
        print(f'Converting a copy of hvgX matrix to {working_dtype_float}...')
        hvg_cp = adata_qc.obsm['hvgX'].copy().astype(working_dtype_float)
        adata_qc.obsm['scaledX'], means, stds = scale_array__dtype_mod(hvg_cp,
                zero_center=True,
                max_value=hvg_scaling_cap,
                copy=True,
                return_mean_std=True,
                working_dtype_float=working_dtype_float)
        # 04b.xiv. Convert scaled count matrix to numpy array so it can be written
        # to hd5 format
        adata_qc.obsm['scaledX'] = np.asarray(adata_qc.obsm['scaledX'])
        # 04b.xv. Because we know that the scaled counts should be much closer
        # to zero, we reduce this matrix to float16 to reduce memory loading.
        adata_qc.uns['hvg_scale_means'] = means.copy()
        adata_qc.uns['hvg_scale_stds'] = stds.copy()
        del means, stds
        adata_qc.obsm['X_pca'] = sc.tl.pca(adata_qc.obsm['scaledX'],
                random_state=numpy_random_state)
        sc.pp.neighbors(adata_qc, 
                        use_rep='X_pca', 
                        n_pcs=n_pcs,
                        method='umap', 
                        metric='cosine',
                        n_neighbors=n_neighbors,
                        random_state=numpy_random_state)
        # 04b.xvi. Perform knn-clustering-based manifold detection and use this
        # to perform nonlinear dimensionality reduction using UMAP, and
        # also Leiden clustering
        sc.tl.umap(adata_qc,
                random_state=numpy_random_state)
        # AMW, 2023-03-09:
        # Note: the current implementation of scanpy.tl.leiden requires
        # a random *seed*; it cannot accept an actual numpy.random.RandomState
        # instance.
        sc.tl.leiden(adata_qc,
                resolution=leiden_resolution,
                random_state=get_rand_seed(numpy_random_state))
        return adata_qc

# 05. Specify hyperparameters to be used for QC and other preprocessing

# 05a. Set thresholds for QC 
# (Note: these were set by following guidelines for removal thresholds
# set out in both CellBender and Seurat documentation, and by inspecting the distributions of these
# three metrics in a subsampled population of our aggregated expression data to adjust these
# guidelines such that cutoffs were meaningful in our single-nucleus data).
# These were also set considering previous observations about human postmortem
# substantia nigra data, noting that some neuron types have low numbers of genes
# per nucleus
# Subramanian et al. 2022 https://pubmed.ncbi.nlm.nih.gov/36575523/
# Welch et al. 2019 https://www.sciencedirect.com/science/article/pii/S0092867419305045
MITO_FRAC_UPPER_PCT = qc_hyperparam_dict['mito_frac_upper_pct']
if not FILTER_ON_MITO_READS:
    MITO_FRAC_UPPER_PCT = 100
print(f'Upper threshold on filtering by mitochondrial read fraction: {MITO_FRAC_UPPER_PCT}')
U_GENE_UPPER_PCT = qc_hyperparam_dict['u_gene_upper_pct']
U_GENE_LOWER_N = qc_hyperparam_dict['u_gene_lower_pct']
N_MOL_UPPER_PCT = qc_hyperparam_dict['n_mol_upper_pct']

# 05b. Also set the number of PCs to use in UMAP calculations
N_PCS = qc_hyperparam_dict['n_pcs']
LEIDEN_RES = qc_hyperparam_dict['leiden_res']

print(f'hyperparameters:\nMITO_FRAC_UPPER_PCT: {MITO_FRAC_UPPER_PCT}\nU_GENE_UPPER_PCT: {U_GENE_UPPER_PCT}\n' + \
        f'U_GENE_LOWER_N: {U_GENE_LOWER_N}\nN_MOL_UPPER_PCT: {N_MOL_UPPER_PCT}\nN_PCS: {N_PCS}\nLEIDEN_RES: {LEIDEN_RES}\n')

# 05c. Set the precision of the accumulators used to compute
# means and variances for scaling such that these do not
# overflow when we compute these values for HVG matrices
# with memory-saving dtypes like uint16/float16.
ACCUMULATOR_DTYPE = 'float32'
if int(working_precision) > 32:
    ACCUMULATOR_DTYPE = 'float64'

# 06. Perform QC on the aggregated data in full or
# pass-fail-only mode, depending on what was specified
# (Note: 25 PCs was determined to be sufficient to capture
# pretty much all of the variance in our hvg expression data,
# based on inspections done in the subsampled data.)
if PASS_FAIL_QC_ONLY:
    print(f'Running QC filtration pass-fail designation only for expression data.')
    adatas_human_qc_pre_harmony_pf_only = qc_and_typical_workflow(adatas_human.copy(),
            umi_percentile=N_MOL_UPPER_PCT,
            gene_percentile=U_GENE_UPPER_PCT,
            n_gene_min=U_GENE_LOWER_N,
            mito_frac_percentile=MITO_FRAC_UPPER_PCT,
            hvg_scaling_cap=MAX_VAL_HVG_SCALING,
            n_hvgs_to_keep=N_HVGS,
            n_pcs=N_PCS,
            n_neighbors=N_NEIGHBORS,
            leiden_resolution=LEIDEN_RES,
            working_dtype_int=working_dtype_int,
            working_dtype_float=working_dtype_float,
            scaling_accumulator_dtype=ACCUMULATOR_DTYPE,
            pass_fail_only_mode=PASS_FAIL_QC_ONLY)
    ahc_pf_only_out_fn = f'{preprocessed_data_dir}/{output_data_pf_only_fn}'
    adatas_human_qc_pre_harmony_pf_only.write(ahc_pf_only_out_fn)
    sys.exit()
else:
    adatas_human_qc_pre_harmony = qc_and_typical_workflow(adatas_human.copy(),
            umi_percentile=N_MOL_UPPER_PCT, 
            gene_percentile=U_GENE_UPPER_PCT, 
            n_gene_min=U_GENE_LOWER_N,
            mito_frac_percentile=MITO_FRAC_UPPER_PCT,
            hvg_scaling_cap=MAX_VAL_HVG_SCALING,
            n_hvgs_to_keep=N_HVGS,
            n_pcs=N_PCS,
            n_neighbors=N_NEIGHBORS,
            leiden_resolution=LEIDEN_RES,
            working_dtype_int=working_dtype_int,
            working_dtype_float=working_dtype_float,
            scaling_accumulator_dtype=ACCUMULATOR_DTYPE)

    del adatas_human

    # 07. Run Harmony for removal of batch effects across datasets.
    # Here, we will use the pool ID as the variable of interest,
    # to capture batch effects without removing true inter-patient
    # variation
    # Harmony computes corrected PCA representations of the count
    # data; see the development .ipynb file for this script
    # (/sc/arion/projects/motor/WILSOA28/demuxlet_analysis/
    # demuxlet__batch_1/terra/cellbender_output_qc__poolwise.ipynb)
    # for more detailed information.
    harmony_output = hm.run_harmony(data_mat=adatas_human_qc_pre_harmony.obsm['X_pca'].copy(),
                                    meta_data=adatas_human_qc_pre_harmony.obs.copy(),
                                    vars_use=harmony_batch_var,
                                    max_iter_harmony=20,
                                    random_state=get_rand_seed(numpy_random_state)
                                    )
    # 07a. Back-calculated scaled (and truncated) expression values in the gene expression
    # basis (for HVGs only) from the Harmonized PCA matrix (in which the coefficients
    # of PCs are adjusted to remove batch effects but the PCs themselves don't change),
    # and use this to construct a new anndata object.
    adatas_human_qc_pc_info = sc.tl.pca(adatas_human_qc_pre_harmony.obsm['scaledX'].copy(),
            dtype='float16',
            return_info=True,
            random_state=numpy_random_state)
    pc_loadings = adatas_human_qc_pc_info[1].copy().astype('float16')
    del adatas_human_qc_pc_info
    harmony_adjusted_X_pca = harmony_output.Z_corr.copy().transpose()
    harmony_adjusted_scaledX_hvg = np.matmul(
            harmony_adjusted_X_pca.copy(),
            pc_loadings)
    del pc_loadings
    # 07b. Undo scaling transforms from hvg data now that it has been corrected
    # in scaled form
    # Do this rescaling at half precision to save space (i.e. when using a large
    # number of HVGs)
    means_hvg_scale = adatas_human_qc_pre_harmony.uns['hvg_scale_means'].copy()
    stds_hvg_scale = adatas_human_qc_pre_harmony.uns['hvg_scale_stds'].copy()
    harmony_adjusted_X_hvg = np.add(
            np.multiply(
                harmony_adjusted_scaledX_hvg.astype('float16'),
                stds_hvg_scale.astype('float16')
                ).astype('float16'),
            means_hvg_scale.astype('float16')) #.astype(working_dtype_float)

    del means_hvg_scale, stds_hvg_scale

    # 07c. Construct a new anndata object with Harmony-corrected hvg expression
    # values, and use this data to compute batch-corrected nearest-neighbor clustering
    # manifold, UMAP embedding, and Leiden clustering (using the same hyperparameters
    # as used in the workflow for the pre-Harmony data).
    adatas_human_qc = ad.AnnData(X=harmony_adjusted_X_hvg.copy(),
                                 dtype=adatas_human_qc_pre_harmony.X.dtype)
    del harmony_adjusted_X_hvg
    adatas_human_qc.obs = adatas_human_qc_pre_harmony.obs.copy()
    if 'leiden' in adatas_human_qc.obs.columns:
        del adatas_human_qc.obs['leiden']
    adatas_human_qc.var = adatas_human_qc_pre_harmony.var.loc[
            adatas_human_qc_pre_harmony.var['highly_variable']
            ].copy()
    adatas_human_qc.uns['n_pass_nuclei_total'] = adatas_human_qc_pre_harmony.uns['n_pass_nuclei_total'].copy()
    adatas_human_qc.uns['hvg'] = adatas_human_qc_pre_harmony.uns['hvg'].copy()
    adatas_human_qc.obsm['X_pca'] = harmony_adjusted_X_pca.copy()
    adatas_human_qc.obsm['scaledX'] = harmony_adjusted_scaledX_hvg.copy()
    sc.pp.neighbors(adatas_human_qc,
                    use_rep='X_pca', 
                    n_pcs=N_PCS,
                    method='umap', 
                    metric='cosine',
                    n_neighbors=N_NEIGHBORS,
                    random_state=numpy_random_state)
    sc.tl.umap(adatas_human_qc,
            random_state=numpy_random_state)
    sc.tl.leiden(adatas_human_qc,
            resolution=LEIDEN_RES,
            random_state=get_rand_seed(numpy_random_state))

    del harmony_adjusted_X_pca, harmony_adjusted_scaledX_hvg

    # 08. Generate 2D UMAP plots of Leiden clustering for both datasets
    # 08a. Pre-Harmony
    plt.figure()
    ax = plt.subplot(1,1,1)
    sc.pl.embedding(adatas_human_qc_pre_harmony,
                    basis='umap',
                    color='leiden',
                    show=False,
                    ax=ax)
    plt.savefig(
            f'{preprocessed_data_dir}/aggregated_leiden_clusters__pre_Harmony.png',
            dpi=300,
            bbox_inches='tight'
            )
    plt.close('all')
    # 08b. Post-Harmony
    plt.figure()
    ax = plt.subplot(1,1,1)
    sc.pl.embedding(adatas_human_qc,
                    basis='umap',
                    color='leiden',
                    show=False,
                    ax=ax)
    plt.savefig(
            f'{preprocessed_data_dir}/aggregated_leiden_clusters.png',
            dpi=300,
            bbox_inches='tight'
            )
    plt.close('all')

    # 09. Save the cleaned, QCed, aggregated,
    # anndata objects to file for downstream cluster typing and differential expression analysis
    # 09a. Set the 'random state' objects to a state object using get_state()
    # or else delete them. These are stored in 
    # uns['neighbors']['params']['random_state']
    # and uns['umap']['params']['random_state']
    # 09a.i. Save pre-Harmony random state params as dictionaries, not tuples
    adatas_human_qc_pre_harmony.uns['neighbors']['params']['random_state'] = adatas_human_qc_pre_harmony.uns[
            'neighbors'][
                    'params'][
                            'random_state'].get_state(
                                    legacy=False
                                    ).copy()
    adatas_human_qc_pre_harmony.uns['umap']['params']['random_state'] = adatas_human_qc_pre_harmony.uns[
            'umap'][
                    'params'][
                            'random_state'].get_state(
                                    legacy=False
                                    ).copy()
    # 09a.ii. Convert post-Harmony random states from tuples to dictionaries
    adatas_human_qc.uns['neighbors']['params']['random_state'] = adatas_human_qc.uns[
            'neighbors'][
                    'params'][
                            'random_state'].get_state(
                                    legacy=False
                                    ).copy()
    adatas_human_qc.uns['umap']['params']['random_state'] = adatas_human_qc.uns[
            'umap'][
                    'params'][
                            'random_state'].get_state(
                                    legacy=False
                                    ).copy()
    # 09b. Save the QCed datasets to file.
    # 09b.i. Pre-Harmony
    ah_pre_harmony_out_fn = f'{preprocessed_data_dir}/{output_data_no_batch_corr_fn}'
    adatas_human_qc_pre_harmony.write(ah_pre_harmony_out_fn)
    # 09b.ii. Post-Harmony
    ah_out_fn = f'{preprocessed_data_dir}/{output_data_batch_corr_fn}'
    adatas_human_qc.write(ah_out_fn)

    # 10. Save the random number generator state to a separate file.
    numpy_random_state_dict = numpy_random_state.get_state(
            legacy=False
            )
    numpy_random_state_output_fn = f'{random_state_dir}/{random_state_fn}'
    with open(numpy_random_state_output_fn,'wb') as outfile:
        pickle.dump(numpy_random_state_dict,
                outfile)

    # 11. Check cluster coverage by pool and patient to make sure no clusters
    # are dominated by a single pool or patient post-Harmony (which would
    # be indicative of persistent barcode-associated artifacts).
    # Also check cluster coverage per pool and patient.
    def check_cluster_coverage(adata_qc: ad.AnnData):
        # 11a. Find the number of pools contributing to each cluster
        cluster_info = pd.DataFrame(columns=['Cluster','N_Patients','N_Pools','Potential_Artifact_1y0n'])
        ci_idx = 0
        for key,value in adata_qc.obs.groupby('leiden'):
            n_pts = len(np.unique(value['patient_ID'].values))
            n_pools = len(np.unique(value['pool_name'].values))
            cluster_info.at[ci_idx,'Cluster'] = key
            cluster_info.at[ci_idx,'N_Patients'] = n_pts
            cluster_info.at[ci_idx,'N_Pools'] = n_pools
            if ((n_pts == 1) or (n_pools == 1)):
                cluster_info.at[ci_idx,'Potential_Artifact_1y0n'] = 1
            else:
                cluster_info.at[ci_idx,'Potential_Artifact_1y0n'] = 0
            ci_idx += 1
        # 11b. Find the number of potentially valid clusters individual patients populate
        potentially_valid_clusters = cluster_info.loc[cluster_info['Potential_Artifact_1y0n']==0]['Cluster'].values.tolist()
        n_pvcs = len(potentially_valid_clusters)
        n_all_clusters = len(cluster_info)
        print(f'Of all {n_all_clusters} clusters, {n_pvcs} are constituted by patients from >1 pool and are considered potentially valid.')
        patient_cluster_info = pd.DataFrame(columns=['Patient_ID','N_Valid_Clusters_Populated','Potential_Issue_1y0n'])
        pci_idx = 0
        for key,value in adata_qc[
                adata_qc.obs['leiden'].isin(potentially_valid_clusters),:
                ].obs.groupby('patient_ID'):
            n_valid_clusters = len(np.unique(value['leiden'].values))
            patient_cluster_info.at[pci_idx,'Patient_ID'] = key
            patient_cluster_info.at[pci_idx,'N_Valid_Clusters_Populated'] = n_valid_clusters
            pci_idx += 1
        # 11c. First, we plot the distribution of number of clusters occupied for the full dataset:
        plt.figure()
        binmax = np.max(patient_cluster_info['N_Valid_Clusters_Populated'].values.tolist())+2
        binedges = [_-0.5 for _ in np.arange(0,binmax,1).tolist()]
        ns_hist, bins_hist, patches = plt.hist(patient_cluster_info['N_Valid_Clusters_Populated'].values.tolist(),
                bins=binedges,
                alpha=0.7,
                color='gray')
        y_max = np.max(ns_hist)+5
        # 11d. Make y_max an even value
        y_max += y_max%2
        pctiles = np.arange(5,100,5)
        cl_pct = np.percentile(patient_cluster_info['N_Valid_Clusters_Populated'].values.tolist(),
                           pctiles)
        for i in range(len(pctiles)):
            plt.plot([cl_pct[i],cl_pct[i]], [0,25], color='crimson', alpha=0.5)
            text_ys = [18]*len(cl_pct)
            # 11d.i. Alternate y location to avoid overlap
            text_y_loc = [18,22][i%2]
            plt.text(x=cl_pct[i]+0.04, y=text_y_loc, s=f'{pctiles[i]}', fontsize='x-small', rotation=45.0)
        plt.xlabel('Number of Valid Clusters Populated by Patient')
        plt.ylabel('Frequency')
        plt.xticks(np.arange(0,binmax)+1,rotation=90)
        plt.yticks(np.arange(0,y_max,2))
        plt.savefig(
                f'{preprocessed_data_dir}/dist_n_valid_clusters_populated_per_pt.png',
                dpi=300,
                bbox_inches='tight'
                )
        plt.close('all')
        # 11e. Flag patients with fewer than some threshold of cluster occupation as potentially problematic.
        threshold_pctile = 5
        threshold_value = np.ceil(cl_pct[np.argwhere(pctiles==threshold_pctile).squeeze().tolist()])
        for idx in patient_cluster_info.index:
            n_clusters = patient_cluster_info.loc[idx]['N_Valid_Clusters_Populated']
            if n_clusters <= threshold_value:
                patient_cluster_info.at[idx,'Potential_Issue_1y0n'] = 1
            else:
                patient_cluster_info.at[idx,'Potential_Issue_1y0n'] = 0
        return cluster_info, patient_cluster_info 

    # 11f. Run this clustering on the aggregated human data
    ah_cluster_info, ah_pt_cluster_info = check_cluster_coverage(adatas_human_qc)
    # 11g. Save cluster population information for inspection
    ah_cluster_info.to_csv(f'{preprocessed_data_dir}/human_aggregated_leiden_cluster_composition.csv',
            index=False)
    ah_pt_cluster_info.to_csv(f'{preprocessed_data_dir}/human_aggregated_pt_population_valid_clusters.csv',
            index=False)



