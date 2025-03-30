import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import argparse
import os
import sys
import yaml
from cellbender.remove_background.downstream import anndata_from_h5 

# This script aggregates a specified list of CellRanger-aligned
# single-nucleus RNA sequencing (snRNA-seq) datasets.
# It also applies filters to each dataset to retain only
# CellBender-cleaned barcodes that are assigned by demuxlets as
# singlets (i.e. as belonging to a single individual), and
# saves the resulting otherwise-raw dataset to file
# for further preprocessing.

# Script setup

# 00. Create argparse object for reading input arguments
parser = argparse.ArgumentParser(description='Filtering and aggregation of expression data.')
parser.add_argument(f'--config-yaml-path',type=str,required=True)
args = parser.parse_args()

# 01. Get analysis parameters from configuration file
cfg = None
with open(args.config_yaml_path) as f:
    cfg = yaml.safe_load(f)
# 01a. Dataset information
# 01a.i. Data parent directory
qc_root_dir = cfg.get('root_dir')
# 01a.ii. Directory with CellBender results per
# dataset (pooled-sample sequencing library)
dataset_CellBender_dir = qc_root_dir + f'/' + cfg.get('dataset_CB_out_dir')
# 01a.iii. Directory with demuxlet singlet files per dataset
dataset_demux_pass_dir = qc_root_dir + f'/' + cfg.get('dataset_demux_out_dir')
# 01a.iv. Directory with demuxlet barcode type
# (singlet, doublet, and ambiguous) counts per dataset
dataset_demux_sng_dbl_amb_dir = qc_root_dir + f'/' + cfg.get('dataset_demux_sng_dbl_amb_counts_dir')
# 01a.v. Dataset output directory
qc_output_dir = qc_root_dir + f'/' + cfg.get('preprocessed_data_dir')
if not os.path.isdir(qc_output_dir):
    os.system(f'mkdir -p {qc_output_dir}')
# 01a.vi. Dataset output file name
qc_output_fn = cfg.get('aggregated_raw_expr_data_fn')
# 01a.vii. Patients to remove from downstream processing and analysis
pt_IDs_to_remove = cfg.get('patients_to_remove')
# 01a.viii. Mapping to reformat pool names used in demuxlet
# so they are consistent with expression data pool names
demux_pool_name_remap = cfg.get('sample_pool_demux_rename_map')

# 02. Read CellBender-filtered count data by pool
output_files = [
        f'{dataset_CellBender_dir}/{_}'
        for _ in
        os.listdir(dataset_CellBender_dir)
        ]
adatas = {}
for out_f in output_files:
    # 02a. Check whether files exist before proceeding
    if os.path.isfile(out_f):
        poolname = out_f.split('/')[-1].split('_')[0]
        extension = out_f.split('.')[-1]
        # 02b. Load expression data files
        adata_curr = None
        if extension=='h5ad':
            adata_curr = ad.read_h5ad(out_f)
        elif extension=='h5':
            adata_curr = anndata_from_h5(file=out_f)
        # 02c. Make sure gene names are unique
        adata_curr.var_names_make_unique() 
        adatas[f'{poolname}']=adata_curr.copy()
        del adata_curr

# 03. Read in demuxlet-passing singlet barcodes by pool
demux_pass_barcode_full_fns = [
        f'{dataset_demux_pass_dir}/{_}'
        for _ in
        os.listdir(dataset_demux_pass_dir)
        ]
demuxlet_barcodes = {}
for full_dpb_fn in demux_pass_barcode_full_fns:
    if os.path.isfile(full_dpb_fn):
        poolname = full_dpb_fn.split('.')[0].split('_')[-1]
        dpb_df = pd.read_csv(full_dpb_fn,
                             index_col=0,
                             dtype={
                                 'BARCODE':str,
                                 'SNG.BEST.GUESS.REFORMATTED':str
                                    }
                            )
        demuxlet_barcodes[poolname] = dpb_df
        del dpb_df

# 04. Pull a list of the (multi-individual) doublets
# identified by demuxlet for each pool, to help 
# estimate the expected number of single-individual
# doublets for Scrublet
demux_n_doublets_all = pd.DataFrame()
demux_barcode_type_counts_full_fns = [
        f'{dataset_demux_sng_dbl_amb_dir}/{_}'
        for _ in
        os.listdir(dataset_demux_sng_dbl_amb_dir)
        ]
for demux_barcode_type_counts_full_fn in demux_barcode_type_counts_full_fns:
    demux_bt_counts_df = pd.read_csv(
            demux_barcode_type_counts_full_fn,
            index_col=0
            )
    if len(demux_n_doublets_all)==0:
        demux_n_doublets_all = demux_bt_counts_df[['Pool_ID','N_DBL_expected_patients']].copy()
    else:
        demux_n_doublets_all = pd.concat(
                [
                    demux_n_doublets_all,
                    demux_bt_counts_df[['Pool_ID','N_DBL_expected_patients']].copy()
                    ],
                ignore_index=True
                )

# 04a. Reformat pool names to be consistent with expression data as needed
for name_old, name_new in demux_pool_name_remap.items():
    demux_n_doublets_all['Pool_ID'] = [
            (name_new + _.split(name_old)[-1])
            if name_old in _
            else _
            for _ in demux_n_doublets_all['Pool_ID'].values.tolist()
            ]

# 04b. Save multi-individual doublet counts to file
demux_dblt_counts_out_fn = f'{qc_output_dir}/demuxlet_multi_ind_doublet_counts.csv'
demux_n_doublets_all.to_csv(demux_dblt_counts_out_fn,
        index=True)

# 05. Filter AnnData objects to keep only demuxlet singlet barcodes
# 05a. Specify a helper function for filtering
def filter_anndata(adatas,inplace_subsample=False,ss_fraction=0.10):
    # 05a.i. Dictionary for filtered, pooled data
    adatas_filt_ppool = {}
    for adkey in adatas.keys():
        adata_curr = adatas[adkey]
        demux_curr = demuxlet_barcodes[adkey]
        # 05a.ii. Record total gene counts for the current pool
        adata_curr.var['n_counts_all_cells'] = np.array(adata_curr.X.sum(axis=0), dtype=int).squeeze()
        # 05a.iii. Filter expression data to keep demuxlet-passing singlets
        adata_pool = adata_curr[
            adata_curr.obs_names.isin(
                demux_curr['BARCODE'].values.tolist()
            )
        ].copy()
        adata_pool.var['n_counts_all_cells_filtered'] = np.array(adata_pool.X.sum(axis=0), dtype=int).squeeze()
        # 05a.iv. Store the total number of demuxlet-passing barcodes per patient as a dictionary
        pool_counts_dict = {}
        demux_pt_groups = demux_curr.groupby(['SNG.BEST.GUESS.REFORMATTED'])['BARCODE'].groups
        for key in demux_pt_groups.keys():
            # 05a.v. Get demuxlet-passing barcodes for the current donor
            pass_barcodes = demux_curr.loc[
                demux_pt_groups[key]
            ]['BARCODE'].values.tolist()
            # 05a.vi. Add the demuxlet-assigned donor ID as an annotation to those barcodes
            # in the filtered expression data
            adata_pt_obs = adata_pool[adata_pool.obs_names.isin(pass_barcodes)].obs_names.values.tolist()
            adata_pool.obs.loc[adata_pt_obs,'patient_ID'] = key
            # 05a.vii. Store the number of demuxlet-passing barcodes for the current donor and pool
            pool_counts_dict[key] = adata_curr[adata_curr.obs_names.isin(pass_barcodes)].shape[0]

        # 05a.viii. Add dictionary with total number of barcodes per patient to the
        # unstructured annotations for the current pool's expression data
        adata_pool.uns['n_pass_nuclei_total'] = pool_counts_dict
        # 05a.ix. If specified, generate a subsampled version of the expression data to return
        if inplace_subsample:
            sc.pp.subsample(adata_pool,
                    fraction=ss_fraction,
                    copy=False,
                    random_state=numpy_random_state)
        # 05a.x. Add expression data for the current pool to the pool anndata dictionary
        adatas_filt_ppool[adkey] = adata_pool
    return adatas_filt_ppool

# 05b. Run filtering on input datasets
adatas_filt_ppool = filter_anndata(adatas, inplace_subsample=False)

# 05c. Pull only data mapped to the human genome from each expression data object (one per pool)
def create_combined_adata_dict(collection_list,inspect_genomes=0,human_only=0):
    adatas_all = {}
    for collection in collection_list:
        for key in collection.keys():
            if len(collection[key]) > 0:
                adata = ad.AnnData(X=collection[key].X.copy(), dtype=collection[key].X.dtype)
                adata.obs = collection[key].obs[['patient_ID']].copy()
                adata.obs['pool_name'] = [key]*adata.X.shape[0]
                adata.var = collection[key].var[
                        ['ambient_expression','feature_type','genome','gene_id','n_counts_all_cells','n_counts_all_cells_filtered']
                        ].copy()
                adata.uns['n_pass_nuclei_total'] = collection[key].uns['n_pass_nuclei_total'].copy()
                if inspect_genomes:
                    # 05c.i. Check total expression for the genomes with material present in each pool
                    # (this step is meant to identify any reads from other genomes that may be contaminating
                    # barcodes assigned to individual donors)
                    genomes = np.unique(adata.var['genome'].values).tolist()
                    print(f'genomes in pool {key}: {genomes}')
                    for genome in genomes:
                        # 05c.ii. Get the expression levels for genes mapped to the current genome
                        exp_genome = adata[:,adata.var['genome'] == genome].X.sum()
                        print(f'expression for {genome}: {exp_genome}')
                    print('\n')
                if human_only:
                    # 05c.ii. Filter the dataset to include only reads from the human genome
                    adata = adata[:,adata.var['genome'] == 'GRCh38'].copy()
                adatas_all[key] = adata
                del adata
            else:
                print(f'Pool {key} has no data; this pool will not be included in downstream QC or analysis.')
    return adatas_all

dataset_list = [adatas_filt_ppool]

adatas_all_human = create_combined_adata_dict(dataset_list,
        inspect_genomes=0,
        human_only=1)

# 05d. Delete individual raw and filtered anndata dictionaries to free up space
del adatas, adatas_filt_ppool

# 05e. Concatenate the filtered, human-only anndata objects
adatas_human = ad.concat(adatas_all_human,join='outer',index_unique='_',fill_value=0,merge='same',uns_merge=None)
# 05e.i. Explicitly build a combined 'n_pass_nuclei_total' dictionary so that the total number of nuclei per patient is preserved
n_pass_nuclei_total_all = {}
for pool_key_token, adata_token in adatas_all_human.items():
    dict_to_add = adata_token.uns['n_pass_nuclei_total']
    for pool_pt_key, n_nuc_pt_pool in dict_to_add.items():
        if pool_pt_key in n_pass_nuclei_total_all.keys():
            n_pass_nuclei_total_all[pool_pt_key] += n_nuc_pt_pool
        else:
            n_pass_nuclei_total_all[pool_pt_key] = n_nuc_pt_pool
adatas_human.uns['n_pass_nuclei_total'] = n_pass_nuclei_total_all.copy()

# 05f. Generate the total number of counts across all cells and in all pools for the concatenated, human-only dataset
def merge_total_gene_counts(adatas, count_col):
    ad_keys = list(adatas.keys())
    gene_counts = pd.DataFrame(data={f'{count_col}_{ad_keys[0]}':adatas[ad_keys[0]].var[count_col].copy()})
    for ad_key in ad_keys[1:]:
        df_next = pd.DataFrame(data={f'{count_col}_{ad_key}':adatas[ad_key].var[count_col].copy()})
        gene_counts = gene_counts.join(df_next)
    gene_counts[count_col] = gene_counts.sum(axis=1)
    return gene_counts[[count_col]]

total_gene_counts_merged = merge_total_gene_counts(adatas_all_human,'n_counts_all_cells')
total_filtered_gene_counts_merged = merge_total_gene_counts(adatas_all_human,'n_counts_all_cells_filtered')
adatas_human.var = adatas_human.var.join(total_gene_counts_merged)
adatas_human.var = adatas_human.var.join(total_filtered_gene_counts_merged)

# 05g. Remove any patients with IDs in the list of patients to remove.
pt_IDs_to_keep = [_ for _ in
        list(
            np.unique(
                adatas_human.obs['patient_ID'].values.tolist()
                )
            )
        if _ not in pt_IDs_to_remove]
adatas_human = adatas_human[adatas_human.obs['patient_ID'].isin(pt_IDs_to_keep)].copy()

# 05h. Set the datatype of adatas_human to one that occupies less memory to save space
# NOTE: the maximum value that can be stored in a uint 16 is ~65,000, so check that
# the max value does not exceed this. If it does, store the count matrix as dtype uint32.
max_count_val_for_dtype = np.max(adatas_human.X.copy())
working_precision =''
if max_count_val_for_dtype < np.iinfo(np.uint32).max:
    working_precision = '32'
else:
    print(f'This count matrix requires 64-bit precision, and many preprocessing steps',
            f' (such as HVG detection) may not function as a result.')
    working_precision = '64'
working_dtype_int = f'uint{working_precision}'
working_dtype_float = f'float{working_precision}'
print(f'Integer dtype for processing: {working_dtype_int}\nFloat dtype for processing: {working_dtype_float}')
adatas_human.X = adatas_human.X.copy().astype(working_dtype_int)

# 05i. Delete human anndata pieces used to created the merged dataset
del adatas_all_human

# 06. Save the demuxlet-/CellBender- filtered, combined, but
# otherwise unprocessed dataset to file.
ahc_filt_comb_out_fn = f'{qc_output_dir}/{qc_output_fn}'
adatas_human.write(ahc_filt_comb_out_fn)

print('\n')

sys.exit()

