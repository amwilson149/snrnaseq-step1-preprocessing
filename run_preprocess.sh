#!/bin/bash
#BSUB -P acc_motor # project name

set -ev

# 01. Set up environment
exec_dir=$( pwd )
cd "${exec_dir}"
prpr_dir="${exec_dir}/scripts"

# 02. Set up config files
# 02a. Specify config file path
cfg="${exec_dir}/configs/config_data_build.yaml"
echo "${cfg}"
# 02b. Add root directory to config
# file if it does not exist in there
# yet (meaning the script hasn't been
# run before)
if ! $( grep -q 'root_dir' ${cfg} ); then
	echo "Initializing config file with current directory as root directory"
	echo "root_dir: '${exec_dir}'" >> ${cfg}
else
	echo "Config file already contains root directory; proceeding"
fi

# 03. Run preprocessing
# 03a. Aggregate CellBender- and demuxlet-processed sequencing data
echo -e "Aggregating datasets...."
python "${prpr_dir}/01a_cellbender_output_qc__build_anndata.py" --config-yaml-path ${cfg}
# 03b. Use Scrublet to assign (single-individual) doublet scores to barcodes
echo -e "Scoring for single-individual doublets...."
python "${prpr_dir}/01b_cellbender_output_qc__find_si_doublets.py" --config-yaml-path ${cfg}

# Stop here and manually validate Scrublet thresholds following Scrublet documentation, then
# proceed with scripts

# 03c. Remove single-individual doublets using validated Scrublet thresholds
echo -e "Removing single-individual doublets using manually validated Scrublet thresholds...."
python "${prpr_dir}/01c_cellbender_output_qc__remove_si_doublets.py" --config-yaml-path ${cfg}
# 03d. Perform QC and preprocessing
echo -e "Performing QC and preprocessing...."
python "${prpr_dir}/01d_cellbender_output_qc__apply_qc_pf_option.py" --config-yaml-path ${cfg}
# 03e. Finalize dataset (ensure it is count-like)
echo -e "Ensuring that data is count-like...."
python "${prpr_dir}/01e_preprocess_qced_batch_corr_count_data.py" --config-yaml-path ${cfg}

# 04. Inspect the number of cell type marker genes among HVGs vs. the number of HVGs retained
# (used to select a number of HVGs to keep for preprocessing)
echo -e "Producing figures inspecting cell type marker gene presence vs. number of HVGs retained"
python "${prpr_dir}/01f_inspect_hvg_var_and_mg_overlap.py" --config-yaml-path ${cfg}


