#!/bin/bash

source /opt/exp_software/neutrino/al9/env.sh
#source /opt/exp_software/neutrino/al9/ROOT/v6.32.06/v6.32.06_install/bin/thisroot.sh (non provato funziona il primo) 

#root -b -q /storage/gpfs_data/neutrino/users/croselli/root_macros/test_snap.cpp (non carica il logon)



root -b -q '/storage/gpfs_data/neutrino/users/croselli/root_macros/rootlogon.C' -q '/storage/gpfs_data/neutrino/users/croselli/root_macros/test_snap.cpp'


