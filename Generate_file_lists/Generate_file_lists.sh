#!/bin/bash

# This script generates file lists for simulation productions in SAND.
# It creates lists of reco, gtrac, edep, and digit ROOT files for a specified range of run numbers.
# Only runs where the digit file exists are included, ensuring complete processing.
# The script is configurable via PROD_TITLE, base_path, base_name, n_start, and n_end variables.
# Output files are named {PROD_TITLE}_{type}_list.txt and contain absolute paths for downstream analysis.
# A list of available productions and their configurations is provided in the comments for reference.

#---------------------------------------ADDITIONAL FILES NOT USED FOR THESIS BUT UP-TO-DATE -----------------------------------------------------

PROD_TITLE="AGG_inner"
base_path="/storage/gpfs_data/neutrino/users/amenga/prod-scripts/production"
base_name="SAND_innervol_gisimple_1M"
n_start=3013
n_end=5012

>${PROD_TITLE}_reco_list.txt
>${PROD_TITLE}_gtrac_list.txt
>${PROD_TITLE}_edep_list.txt 
>${PROD_TITLE}_digit_list.txt


for i in $(seq ${n_start} ${n_end}); do
  base_dir="${base_path}/${base_name}/${base_name}_${i}"

  percent=$(( (i - ${n_start} + 1) * 100 / (${n_end} - ${n_start} + 1) ))
  printf "\rGenerating file lists: %3d%%" "$percent"

  reco="${base_dir}/sand-events.${i}.reco.root"
  gtrac="${base_dir}/sand-events.${i}.gtrac.root"
  edep="${base_dir}/sand-events.${i}.edep.root"
  digit="${base_dir}/sand-events.${i}.digit.root"

  if [ -f "$digit" ]; then
    echo "$edep" >> ${PROD_TITLE}_edep_list.txt
    echo "$reco" >> ${PROD_TITLE}_reco_list.txt
    echo "$gtrac" >> ${PROD_TITLE}_gtrac_list.txt
    echo "$digit" >> ${PROD_TITLE}_digit_list.txt
  fi



done
printf "\nGenerated lists! \n"

#---------------------------------------------------------FINAL FILES USED FOR THESIS-----------------------------------------------------------

#PROD_TITLE="FINAL_inner"
#base_path="/storage/gpfs_data/neutrino/users/amenga/prod-scripts/production/"
#base_name="SAND_innervol_gisimple_1M"
#n_start=1
#n_end=3012



#**********************************************OLDER PRODUCTIONS - NOT USED ANYMORE******************************************************

#-------------------------------------------------------------MINI-PRODUZIONE----------------------------------------------------------

# PROD_TITLE="MINI"
# base_path="/storage/gpfs_data/neutrino/users/croselli/productions/prod_al9/prod_tests/"
# base_name="SAND_innerVOL_10"
# n_start=0
# n_end=2000

#------------------------------------------------------------------SAND-------------------------------------------------------------------------------------
# PROD_TITLE="SAND"
# base_path="/storage/gpfs_data/neutrino/users/croselli/productions/prod_al9/prod_tests/"
# base_name="SAND"
# n_start=0
# n_end=1000

#---------------------------------------------------------------SAND_400_EVT-------------------------------------------------------------------------------------

# PROD_TITLE="tot_EVT"
# base_path="/storage/gpfs_data/neutrino/users/croselli/productions/prod_al9/prod_tests/"
# base_name="SAND_400EVT"
# n_start=0
# n_end=1000

#--------------------------------------------------------------SAND_400_EVT_2-------------------------------------------------------------------------------------

# PROD_TITLE="tot_EVT"
# base_path="/storage/gpfs_data/neutrino/users/croselli/productions/prod_al9/prod_tests/"
# base_name="SAND_400EVT_2"
# n_start=0
# n_end=3500


#----------------------------------------------------------------SAND_AN-------------------------------------------------------------------------------------

# PROD_TITLE="tot_EVT"
# base_path="/storage/gpfs_data/neutrino/users/croselli/productions/prod_al9/prod_tests/"
# base_name="SAND_AN"
# n_start=0
# n_end=1075



#--------------------------------------------------------------GSIMPLE_TEST-------------------------------------------------------------------------------------

# PROD_TITLE="GSIM"
# base_path="/storage/gpfs_data/neutrino/users/croselli/productions/prod_al9/prod_tests/"
# base_name="SAND_innerVOL_gsimple_test"
# n_start=1
# n_end=20

#------------------------------------------------------------INNER_RECO_LIST----------------------------------

# PROD_TITLE="inner"
# base_path="/storage/gpfs_data/neutrino/users/croselli/productions/prod_al9/prod_tests/"
# base_name="SAND_innerVOL_10"
# n_start=0
# n_end=7250

#---------------------------------------------------------------GS_INNER-----------------------------------------

# PROD_TITLE="GS_inner"
# base_path="/storage/gpfs_data/neutrino/users/amenga/prod-scripts/production/"
# base_name="SAND_innervol_gisimple_1M"
# n_start=1
# n_end=1011
