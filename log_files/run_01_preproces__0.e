
# 01. Set up environment
# 01a. Activate conda environment
ml anaconda3/2020.11
module  load 'anaconda3/2020.11'
Shell debugging temporarily silenced: export LMOD_SH_DBG_ON=1 for Lmod's output

The following have been reloaded with a version change:
  1) gcc/14.2.0 => gcc/8.3.0

Shell debugging restarted
ml -python
module  load '-python'
Shell debugging temporarily silenced: export LMOD_SH_DBG_ON=1 for Lmod's output
Shell debugging restarted
source /hpc/packages/minerva-centos7/anaconda3/2020.11/etc/profile.d/conda.sh
export CONDA_EXE='/hpc/packages/minerva-centos7/anaconda3/2020.11/bin/conda'
export _CE_M=''
export _CE_CONDA=''
export CONDA_PYTHON_EXE='/hpc/packages/minerva-centos7/anaconda3/2020.11/bin/python'

# Copyright (C) 2012 Anaconda, Inc
# SPDX-License-Identifier: BSD-3-Clause

__add_sys_prefix_to_path() {
    # In dev-mode CONDA_EXE is python.exe and on Windows
    # it is in a different relative location to condabin.
    if [ -n "${_CE_CONDA}" ] && [ -n "${WINDIR+x}" ]; then
        SYSP=$(\dirname "${CONDA_EXE}")
    else
        SYSP=$(\dirname "${CONDA_EXE}")
        SYSP=$(\dirname "${SYSP}")
    fi

    if [ -n "${WINDIR+x}" ]; then
        PATH="${SYSP}/bin:${PATH}"
        PATH="${SYSP}/Scripts:${PATH}"
        PATH="${SYSP}/Library/bin:${PATH}"
        PATH="${SYSP}/Library/usr/bin:${PATH}"
        PATH="${SYSP}/Library/mingw-w64/bin:${PATH}"
        PATH="${SYSP}:${PATH}"
    else
        PATH="${SYSP}/bin:${PATH}"
    fi
    \export PATH
}

__conda_hashr() {
    if [ -n "${ZSH_VERSION:+x}" ]; then
        \rehash
    elif [ -n "${POSH_VERSION:+x}" ]; then
        :  # pass
    else
        \hash -r
    fi
}

__conda_activate() {
    if [ -n "${CONDA_PS1_BACKUP:+x}" ]; then
        # Handle transition from shell activated with conda <= 4.3 to a subsequent activation
        # after conda updated to >= 4.4. See issue #6173.
        PS1="$CONDA_PS1_BACKUP"
        \unset CONDA_PS1_BACKUP
    fi

    \local cmd="$1"
    shift
    \local ask_conda
    CONDA_INTERNAL_OLDPATH="${PATH}"
    __add_sys_prefix_to_path
    ask_conda="$(PS1="$PS1" "$CONDA_EXE" $_CE_M $_CE_CONDA shell.posix "$cmd" "$@")" || \return $?
    rc=$?
    PATH="${CONDA_INTERNAL_OLDPATH}"
    \eval "$ask_conda"
    if [ $rc != 0 ]; then
        \export PATH
    fi
    __conda_hashr
}

__conda_reactivate() {
    \local ask_conda
    CONDA_INTERNAL_OLDPATH="${PATH}"
    __add_sys_prefix_to_path
    ask_conda="$(PS1="$PS1" "$CONDA_EXE" $_CE_M $_CE_CONDA shell.posix reactivate)" || \return $?
    PATH="${CONDA_INTERNAL_OLDPATH}"
    \eval "$ask_conda"
    __conda_hashr
}

conda() {
    if [ "$#" -lt 1 ]; then
        "$CONDA_EXE" $_CE_M $_CE_CONDA
    else
        \local cmd="$1"
        shift
        case "$cmd" in
            activate|deactivate)
                __conda_activate "$cmd" "$@"
                ;;
            install|update|upgrade|remove|uninstall)
                CONDA_INTERNAL_OLDPATH="${PATH}"
                __add_sys_prefix_to_path
                "$CONDA_EXE" $_CE_M $_CE_CONDA "$cmd" "$@"
                \local t1=$?
                PATH="${CONDA_INTERNAL_OLDPATH}"
                if [ $t1 = 0 ]; then
                    __conda_reactivate
                else
                    return $t1
                fi
                ;;
            *)
                CONDA_INTERNAL_OLDPATH="${PATH}"
                __add_sys_prefix_to_path
                "$CONDA_EXE" $_CE_M $_CE_CONDA "$cmd" "$@"
                \local t1=$?
                PATH="${CONDA_INTERNAL_OLDPATH}"
                return $t1
                ;;
        esac
    fi
}

if [ -z "${CONDA_SHLVL+x}" ]; then
    \export CONDA_SHLVL=0
    # In dev-mode CONDA_EXE is python.exe and on Windows
    # it is in a different relative location to condabin.
    if [ -n "${_CE_CONDA+x}" ] && [ -n "${WINDIR+x}" ]; then
        PATH="$(\dirname "$CONDA_EXE")/condabin${PATH:+":${PATH}"}"
    else
        PATH="$(\dirname "$(\dirname "$CONDA_EXE")")/condabin${PATH:+":${PATH}"}"
    fi
    \export PATH

    # We're not allowing PS1 to be unbound. It must at least be set.
    # However, we're not exporting it, which can cause problems when starting a second shell
    # via a first shell (i.e. starting zsh from bash).
    if [ -z "${PS1+x}" ]; then
        PS1=
    fi
fi
conda activate cellbender-w-scrublet-env
PS1='(cellbender-w-scrublet-env) '
export PATH='/sc/arion/work/wilsoa28/.conda/envs/cellbender-w-scrublet-env/bin:/hpc/packages/minerva-centos7/anaconda3/2020.11/bin:/hpc/packages/minerva-centos7/gcc/8.3.0_32b/bin:/hpc/packages/minerva-centos7/anaconda3/2020.11/condabin:/hpc/users/wilsoa28/google-cloud-sdk/bin:/hpc/users/wilsoa28/git-filter-repo:/hpc/packages/minerva-rocky9/git/2.46.0/bin:/hpc/packages/minerva-common/vim/8.0/bin:/hpc/lsf/10.1/linux3.10-glibc2.17-x86_64/etc:/hpc/lsf/10.1/linux3.10-glibc2.17-x86_64/bin:/bin:/usr/bin:/usr/mbin:/local/bin:/usr/local:/usr/ucb:/usr/local/sbin:/usr/sbin:/usr/lpp/mmfs/bin:/hpc/users/wilsoa28/.local/bin:/hpc/users/wilsoa28/bin'
export CONDA_PREFIX='/sc/arion/work/wilsoa28/.conda/envs/cellbender-w-scrublet-env'
export CONDA_SHLVL='1'
export CONDA_DEFAULT_ENV='cellbender-w-scrublet-env'
export CONDA_PROMPT_MODIFIER='(cellbender-w-scrublet-env) '
export CONDA_EXE='/hpc/packages/minerva-centos7/anaconda3/2020.11/bin/conda'
export _CE_M=''
export _CE_CONDA=''
export CONDA_PYTHON_EXE='/hpc/packages/minerva-centos7/anaconda3/2020.11/bin/python'
# 01b. Get root directory
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
/sc/arion/work/wilsoa28/.conda/envs/cellbender-w-scrublet-env/lib/python3.8/site-packages/anndata/_core/anndata.py:1830: UserWarning: Variable names are not unique. To make them unique, call `.var_names_make_unique`.
  utils.warn_names_duplicates("var")
/sc/arion/work/wilsoa28/.conda/envs/cellbender-w-scrublet-env/lib/python3.8/site-packages/anndata/_core/anndata.py:1830: UserWarning: Variable names are not unique. To make them unique, call `.var_names_make_unique`.
  utils.warn_names_duplicates("var")
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
Calculating mito_frac, assuming mito genes start with
genome prefix plus "mt-" or "MT-"2415 cells before cell QC
2031 cells after cell QC
Detecting highly variable genes with seurat_v3...
2025-03-30 12:58:32,434 - harmonypy - INFO - Iteration 1 of 20
2025-03-30 12:58:32,618 - harmonypy - INFO - Iteration 2 of 20
2025-03-30 12:58:32,802 - harmonypy - INFO - Iteration 3 of 20
2025-03-30 12:58:32,986 - harmonypy - INFO - Iteration 4 of 20
2025-03-30 12:58:33,170 - harmonypy - INFO - Iteration 5 of 20
2025-03-30 12:58:33,354 - harmonypy - INFO - Iteration 6 of 20
2025-03-30 12:58:33,538 - harmonypy - INFO - Iteration 7 of 20
2025-03-30 12:58:33,723 - harmonypy - INFO - Iteration 8 of 20
2025-03-30 12:58:33,907 - harmonypy - INFO - Iteration 9 of 20
2025-03-30 12:58:34,090 - harmonypy - INFO - Converged after 9 iterations
# 03e. Finalize dataset (ensure it is count-like)
echo -e "Ensuring that data is count-like...."
python "${prpr_dir}/01e_preprocess_qced_batch_corr_count_data.py" --config-yaml-path ${cfg}

# 04. Inspect the number of cell type marker genes among HVGs vs. the number of HVGs retained
# (used to select a number of HVGs to keep for preprocessing)
echo -e "Producing figures inspecting cell type marker gene presence vs. number of HVGs retained"
python "${prpr_dir}/01f_inspect_hvg_var_and_mg_overlap.py" --config-yaml-path ${cfg}



