#!/bin/bash
# This script is designed to run all the necessary scripts in the Solid_Hydrogen2 analysis workflow in batch mode.

# Export paths to run jupyter scripts
export HOME=$PWD
export IPYTHONDIR=$PWD/.ipython
export JUPYTER_CONFIG_DIR=$PWD/.jupyter
export JUPYTER_RUNTIME_DIR=$PWD/.jupyter_runtime
export MPLCONFIGDIR=$PWD/.matplotlib

# Source code environment
source /opt/exp_software/neutrino/al9/env.sh

# Activate your desired python environment 
source /storage/gpfs_data/neutrino/users/battisti/python_venvs/basic_venv_3_9_21/bin/activate

base_path="/storage/gpfs_data/neutrino/users/battisti/hydrogen_analysis_tests/Solid_Hydrogen2"

for folder_name in "Generate_file_lists" "Snap_file_processing" "Jupyter_analysis" "Unfolding"; do
        cd "${base_path}/${folder_name}" || exit 1
        echo "Processing folder: ${base_path}/${folder_name} with files:"
        ls -la

        if [[ -f "${folder_name}.sh" ]]; then
            bash "${folder_name}.sh"
        else
            echo "Error: ${folder_name}.sh not found in ${folder_name}."
            exit 1
        fi

        cd "${base_path}" || exit 1
done
