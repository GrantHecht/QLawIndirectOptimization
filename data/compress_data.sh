#!/bin/bash

# Navigate to base of repository
parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd $parent_path

# First compress data/TAES directory
tar -zcvf ./tar/taes-data.tar.gz ./TAES