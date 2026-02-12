#!/bin/bash

if ! ls *.so 1> /dev/null 2>&1; then
    bash generate_dict.sh
fi

source /opt/exp_software/neutrino/al9/env.sh
root -b -q '/storage/gpfs_data/neutrino/users/battisti/hydrogen_analysis_tests/Solid_Hydrogen2/Snap_file_processing/rootlogon.C' -q '/storage/gpfs_data/neutrino/users/battisti/hydrogen_analysis_tests/Solid_Hydrogen2/Snap_file_processing/Snap_file_processing.cpp'


