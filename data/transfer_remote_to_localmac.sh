#!/bin/sh

# Navigate to base of repository
parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd $parent_path

# Transfer file
scp granthec@traj-optim.eng.buffalo.edu:/data/users/granthec/devel/julia/QLawIndirectOptimization/data/tar/taes-data.tar.gz ./tar

# Unpack into data/TAES
tar -xvzf tar/taes-data.tar.gz 
