import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import scipy.sparse as sp
import matplotlib.pyplot as plt
import os
import sys
import seaborn as sns
import argparse
import yaml

# This script completes snRNA-seq data preprocessing,
# by making adjustments to scaled, batch-corrected,
# dimensionality-reduced data including
# (1) clipping the lower tail of the scaled, adjusted
#     values to remove negative post-batch-correction
#     counts that correspond to zero counts
# (2) shifting gene counts so minimum count value matches
#     the pre-adjustment minimum (typically shifts are <1)
# (3) rescaling the data and clipping low counts at 0
# (4) converting data to integer data type.

# NOTE: The imported, batch-corrected data *has* a naively
# rescaled set of HVG counts as .X; however,
# this rescaling does not account for the above small adjustments
# that would in principle be required to make the data count-like,
# since we have retained dimensionality reduction and
# batch correction adjustments for our analysis.

# This script makes those adjustments, starting from the
# *scaled* count matrix of the batch-corrected expression
# data (in .obsm['scaledX']).
# Non-batch-corrected expression data was used to apply
# these corrections, which is the reason it is used
# as an input to this script.
# Notably, we have found these adjustments to be minor
# but to be important for ensuring, e.g., that count
# data is non-negative.
# For this reason, we recommend that the output from this
# step is used in downstream analyses, rather than the
# output files from the previous step.


# Script setup

# 00. Create argparse object for reading input arguments
parser = argparse.ArgumentParser(description='Rescaling of QCed, preprocessed expression data.')
parser.add_argument(f'--config-yaml-path',type=str,required=True)
args = parser.parse_args()

# 01. Get analysis parameters from configuration file
cfg = None
with open(args.config_yaml_path) as f:
    cfg = yaml.safe_load(f)
# 01a. Dataset information
# 01a.i. Data root directory
qc_root_dir = cfg.get('root_dir')
# 01a.iii. Data directory containing preprocessed expression data files
preprocessed_data_dir = qc_root_dir + f'/' + cfg.get('preprocessed_data_dir')
# 01a.iv. Preprocessed expression data file name
input_data_fn = cfg.get('preproc_exp_data_batch_corr_scaled_fn')
# 01a.v. Preprocessed, NOT-BATCH-CORRECTED, un-rescaled expression
# data, which has scaling values to be used for rescaling of the batch-corrected data.
input_pre_harm_data_fn = cfg.get('preproc_exp_data_no_batch_corr_raw_fn')
# 01a.v. Data output directory for rescaled, preprocessed expression data
qc_preproc_output_data_dir = qc_root_dir + f'/' + cfg.get('post_qc_preproc_data_dir')
if not os.path.isdir(qc_preproc_output_data_dir):
    os.system(f'mkdir -p {qc_preproc_output_data_dir}')
# 01a.vi. Rescaled, preprocessed expression data output file name
qc_output_fn = cfg.get('post_qc_preproc_data_fn')

# 02. Import expression data
# 02a. Batch-corrected, QCed data
ahc_fn = f'{preprocessed_data_dir}/{input_data_fn}'
adatas_human_qc = ad.read_h5ad(ahc_fn)

# 02b. Non-batch-corrected, QCed data
ahc_ph_fn = f'{preprocessed_data_dir}/{input_pre_harm_data_fn}'
adatas_human_qc_ph = ad.read_h5ad(ahc_ph_fn)

# 03. Compute expression data adjustments

# 03a. Define a function to compare adjusted vs. unadjusted
# expression data per HVG
def dist_comparison_test(
        hvg_new,
        hvg_orig,
        hvg_names_raw,
        truncate_data=False,
        pctile_clip=1,
        out_dir='./',
        chunk_size=20,
        max_chunk_index=11
        ):
    # 03a.i. Define a file name tag to indicate
    #whether data being tested was clipped
    fn_clipping_tag=''
    if truncate_data == True:
        fn_clipping_tag = f'_clipped_{pctile_clip}_pctile'
    # 03a.ii. Define a dummy DataFrame for storing data
    # for comparison plots
    df_curr = pd.DataFrame()
    if not os.path.isdir(out_dir):
        os.system(f'mkdir -p {out_dir}')
    chunk_index = 0
    # 03a.iii. Assign dummy chunk index
    chunk_start = 0
    chunk_end = 0
    while chunk_index < max_chunk_index:
        chunk_start = chunk_index*chunk_size
        chunk_end = (chunk_index+1)*chunk_size
        arr_new = hvg_new[:,chunk_start:chunk_end]
        arr_orig = hvg_orig[:,chunk_start:chunk_end]
        hvg_names = [_.split('_______________')[-1] for _ in hvg_names_raw]
        # 03a.iv. Generate original/new joint DataFrame for seaborn plots
        do_exp = []
        do_gene = []
        do_dataset = []
        for c_idx in np.arange(arr_orig.shape[1]):
            # 03a.v. If data truncation is specified,
            # keep only the expression values within the
            # specified percentile range [0+p:100-p]
            col_curr = list(arr_orig[:,c_idx])
            if truncate_data == True:
                pct_lims = np.percentile(col_curr,
                        [pctile_clip,100-pctile_clip])
                col_curr = [_ for _ in col_curr if (_>=pct_lims[0] and _<=pct_lims[1])]
            gene_curr = [hvg_names[c_idx]]*len(col_curr)
            dataset_curr = ['orig']*len(col_curr)
            do_exp.extend(col_curr)
            do_gene.extend(gene_curr)
            do_dataset.extend(dataset_curr)
        df_orig = pd.DataFrame(data={'expression':do_exp,
            'gene_name':do_gene,
            'dataset':do_dataset})
        dn_exp = []
        dn_gene = []
        dn_dataset = []
        for c_idx in np.arange(arr_new.shape[1]):
            col_curr = list(arr_new[:,c_idx])
            if truncate_data == True:
                pct_lims = np.percentile(col_curr,
                        [pctile_clip,100-pctile_clip])
                col_curr = [_ for _ in col_curr if (_>=pct_lims[0] and _<=pct_lims[1])]
            gene_curr = [hvg_names[c_idx]]*len(col_curr)
            dataset_curr = ['new']*len(col_curr)
            dn_exp.extend(col_curr)
            dn_gene.extend(gene_curr)
            dn_dataset.extend(dataset_curr)
        df_new = pd.DataFrame(data={'expression':dn_exp,
            'gene_name':dn_gene,
            'dataset':dn_dataset})
        df_curr = pd.concat([df_orig,df_new],ignore_index=True)
        df_curr['gene_name'] = df_curr['gene_name'].astype('category')
        df_curr['dataset'] = df_curr['dataset'].astype('category')
        df_curr['expression'] = df_curr['expression'].astype('float')
        del df_orig, df_new
        # 03a.vi. Make a split violin plot
        ax = sns.violinplot(data=df_curr,
                x='gene_name',
                y='expression',
                hue='dataset',
                split=True,
                cut=0)
        plt.xticks(rotation=45)
        plt.title(f'HVG Batch {chunk_index}')
        plt.savefig(f'{out_dir}/hvg_scaled_expr_violin{fn_clipping_tag}_{chunk_index:02}.png',
                dpi=300)
        plt.close()
        # 03a.vii. Make a box plot
        ax = sns.catplot(data=df_curr,
                x='gene_name',
                y='expression',
                hue='dataset',
                kind='boxen',
                orient='v')
        plt.xticks(rotation=45)
        plt.title(f'HVG Batch {chunk_index}')
        plt.savefig(f'{out_dir}/hvg_scaled_expr_box{fn_clipping_tag}_{chunk_index:02}.png',
                dpi=300)
        plt.close()
        chunk_index += 1

# 03b. Define a function to optionally run the above comparison test,
# and to perform clipping, shifting, rescaling, and lower-bound zero
# clipping of adjusted, scaled HVG expression data
def preprocess_harmonized_expression_data(adata_bc,
        adata_orig,
        run_test=False,
        truncate_data_test=False,
        pctile_bilateral_clip_test=1,
        pctile_lb_clip=0,
        pctile_ub_clip=100):
    # 03b.i. Get adjusted, scaled HVG expression levels
    # (e.g. adjusted PCA embeddings projected onto PC loadings)
    sc_adj_data_0 = adata_bc.obsm['scaledX'].copy()
    # 03b.ii. Perform test comparisons on original and adjusted
    # expression levels if that option is selected
    hvg_orig = adata_orig.obsm['scaledX'].copy()
    out_dir = f'{qc_preproc_output_data_dir}/hvg_rescaling_tests/'
    hvg_names_raw = adata_orig.var_names
    if run_test == True:
        dist_comparison_test(hvg_new=sc_adj_data_0,
                hvg_orig=hvg_orig,
                hvg_names_raw=hvg_names_raw,
                truncate_data=truncate_data_test,
                pctile_clip=pctile_bilateral_clip_test,
                out_dir=out_dir,
                chunk_size=20,
                max_chunk_index=11)
    else:
        # 03b.iii. Pre-process "raw" adjusted data (i.e. the projection
        # of adjusted HVG PC embeddings onto the gene expression basis)
        # 03b.iv. Clip adjusted expression data at defined upper and lower
        # percentiles (default 0 and 100th percentiles, i.e. no clipping)
        print('Preprocessing HVG expression data...')
        print('Clipping...')
        lb_vals = np.percentile(sc_adj_data_0,
                pctile_lb_clip,
                axis=0)
        ub_vals = np.percentile(sc_adj_data_0,
                pctile_ub_clip,
                axis=0)
        for c_idx in range(sc_adj_data_0.shape[1]):
            # 03b.iii.A. Lower bound clip
            sc_adj_data_0[:,c_idx][
                    sc_adj_data_0[:,c_idx] < lb_vals[c_idx]
                    ] = lb_vals[c_idx]
            # 03b.iii.B. Upper bound clip
            sc_adj_data_0[:,c_idx][
                    sc_adj_data_0[:,c_idx] > ub_vals[c_idx]
                    ] = ub_vals[c_idx]

        # 03b.iv. Shift each adjusted distribution up to the minimum value
        # of the original, scaled HVG distribution (so that no counts
        # are negative after rescaling of expression values)
        print('Shifting...')
        hvg_orig_scaled_minima = np.min(hvg_orig,
                axis=0)
        sc_adj_data_minima = np.min(sc_adj_data_0,
                axis=0)
        deltas = np.subtract(hvg_orig_scaled_minima.copy(),
                sc_adj_data_minima.copy())
        delta_arr = np.tile(deltas,
                (sc_adj_data_0.shape[0],1))
        sc_adj_data_shift = np.add(sc_adj_data_0.copy(),
                delta_arr.copy())
        # 03b.v. Re-scale adjusted expression values using the parameters
        # of the original downscaling
        print('Rescaling...')
        means = adata_orig.uns['hvg_scale_means'].copy()
        stds = adata_orig.uns['hvg_scale_stds'].copy()
        sc_adj_data_rescale = np.add(
                np.multiply(
                    sc_adj_data_shift.astype('float16'),stds.astype('float16')
                    ),
                means.astype('float16'))
        # 03b.v.A. Sanity check: ensure the resulting minimum adjusted
        # expression values are very close to 0, by printing
        # the percentile spread of these minima
        print('Performing sanity check...')
        new_minima_sanity_check = np.min(sc_adj_data_rescale,
                axis=0)
        pctiles_sanity_check = np.arange(0,100,10)
        min_pctile_vals_sanity_check = np.percentile(new_minima_sanity_check,
                pctiles_sanity_check)
        print(f'Percentile values of adjusted HVG minimum expression values:')
        print(f'Percentile\t\tMin. Val.')
        for i,q in enumerate(pctiles_sanity_check):
            print(f'{q}\t\t{min_pctile_vals_sanity_check[i]}')
        # 03b.vi. Clip the resulting expression values such that they all have
        # a hard minimum of 0 (to allow for downstream analysis)
        # Pending sanity check, most of the HVG expression distributions
        # should already have minimum values of 0 by this point anyway
        print('Clipping count values to a zero minimum...')
        sc_adj_data_rescale[
                sc_adj_data_rescale < 0
                ] = 0
        return sc_adj_data_rescale.copy() 

# 04. Run tests to explore how original and adjusted HVG expression
# distributions compare
# 04a. No bilateral clipping of distributions (so we can see all outliers)
#preprocess_harmonized_expression_data(adatas_human_qc,
#        adatas_human_qc_ph,
#        run_test=True)
# 04b. Bilateral clipping of the outer 1st percentiles to explore how
# many outliers this removes and to better view the bulk distributions
#preprocess_harmonized_expression_data(adatas_human_qc,
#    adatas_human_qc_ph,
#    run_test=True,
#    truncate_data_test=True)
# 04c. Bilateral clipping of the outer 5th percentiles
#preprocess_harmonized_expression_data(adatas_human_qc,
#    adatas_human_qc_ph,
#    run_test=True,
#    truncate_data_test=True,
#    pctile_bilateral_clip_test=5)


# 05. Run actual pre-processing of the adjusted, scaled data
# From our explorations, we will do lower bound clipping at
# the 1st percentile and no upper bound clipping
hvgX_adj_preprocessed_float = preprocess_harmonized_expression_data(adatas_human_qc,
        adatas_human_qc_ph,
        pctile_lb_clip=1)
# 05a. Convert float values to integer values to mimic counts (see test below
# for verification that this approach makes sense for our data)
hvgX_adj_preprocessed = hvgX_adj_preprocessed_float.copy().astype('int32')

# 05b. Specify a test to look at the distributions of some
# genes of interest to see what they look like as non-integers vs. integers
def run_integer_conversion_test(adatas_human_qc):
    ahc_part = adatas_human_qc[adatas_human_qc.obs.leiden == '4',:].copy()
    ahcX = ahc_part.X.todense().copy()
    # Look at the distributions of HVGs with among the smallest and largest
    # expression levels in this test group to see how sensitive the distributions
    # might be to these changes in expression value
    gene_means = np.mean(ahcX,
            axis=0).tolist()[0]
    test_idxs = np.argsort(gene_means).tolist()[:5] + np.argsort(gene_means).tolist()[-5:]
    test_out_dir = f'{qc_preproc_output_data_dir}/tests'
    if not os.path.isdir(test_out_dir):
        os.system(f'mkdir -p {test_out_dir}')
    for idx_curr in range(len(test_idxs)):
        test_idx = test_idxs[idx_curr]
        ge_init = ahcX[:,test_idx]
        ge_conv_int = ge_init.astype('int32')
        plt.figure()
        plt.hist(ge_init,color='tab:blue',alpha=0.7,label='rescaled float')
        plt.hist(ge_conv_int,color='tab:red',alpha=0.7,label='rf conv. to int')
        plt.legend()
        plt.savefig(f'{test_out_dir}/test_conversion_to_int__expt_{idx_curr}).png')
        plt.close()

# 05c. Run this test to verify that converting counts to integers
# functions as expected
#adatas_human_qc_test = adatas_human_qc.copy()
#adatas_human_qc_test.X = sp.csr_matrix(hvgX_adj_preprocessed_float.copy())
#run_integer_conversion_test(adatas_human_qc=adatas_human_qc_test)

# 06. Replace input expression count data (.X) with rescaled, adjusted values 
print(f'Shape check:\nOriginal data: {adatas_human_qc.X.shape}\nPre-processed data: {hvgX_adj_preprocessed.shape}')
adatas_human_qc.X = sp.csr_matrix(hvgX_adj_preprocessed.copy())

# 07. Save the preprocessed dataset for downstream analysis
ahc_preprocessed_out_fn = f'{qc_preproc_output_data_dir}/{qc_output_fn}'
adatas_human_qc.write(ahc_preprocessed_out_fn)


