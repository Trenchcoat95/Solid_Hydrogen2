#!/bin/bash

base_path="/storage/gpfs_data/neutrino/users/battisti/hydrogen_analysis_tests/Solid_Hydrogen2/"

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
