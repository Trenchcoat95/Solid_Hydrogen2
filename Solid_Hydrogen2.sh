#!/bin/bash

read -p "Do you want to process the files as well? (y/n) " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    for folder_name in "Generate_file_lists" "Snap_file_processing" "Jupyter_analysis" "Unfolding"; do
            cd "${folder_name}" || exit 1

            if [[ -f "${folder_name}.sh" ]]; then
                bash "${folder_name}.sh"
            else
                echo "Error: ${folder_name}.sh not found in ${folder_name}."
                exit 1
            fi

            cd .. || exit 1
done
else
    echo "Skipping file processing, only running analysis scripts. You can run the individual scripts in each folder later if needed."
    for folder_name in "Jupyter_analysis" "Unfolding"; do
            cd "${folder_name}" || exit 1

            if [[ -f "${folder_name}.sh" ]]; then
                bash "${folder_name}.sh"
            else
                echo "Error: ${folder_name}.sh not found in ${folder_name}."
                exit 1
            fi

            cd .. || exit 1
done
fi