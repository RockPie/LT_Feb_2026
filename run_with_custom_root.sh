#!/bin/bash
# Wrapper script to ensure executables use custom ROOT installation
# This script removes conda from the environment temporarily

# Deactivate conda by removing it from PATH and LD_LIBRARY_PATH
export PATH=$(echo "$PATH" | tr ':' '\n' | grep -v 'anaconda3/envs/general' | tr '\n' ':' | sed 's/:$//')
export LD_LIBRARY_PATH=$(echo "$LD_LIBRARY_PATH" | tr ':' '\n' | grep -v 'anaconda3/envs/general' | tr '\n' ':' | sed 's/:$//')

# Set ROOT paths explicitly  
export ROOTSYS=/home/shihai/sw/root/root_install
export PATH=$ROOTSYS/bin:$PATH
export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH

# Clear conda variables
unset CONDA_PREFIX
unset CONDA_DEFAULT_ENV
unset CONDA_PROMPT_MODIFIER

# Run the requested command
exec "$@"
