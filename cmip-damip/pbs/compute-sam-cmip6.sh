#!/bin/bash

#PBS -P lo70
#PBS -q normal
#PBS -l ncpus=30
#PBS -l mem=32GB
#PBS -l walltime=0:30:00
#PBS -l wd
#PBS -l jobfs=400GB
#PBS -l storage=gdata/ux62+scratch/k10+gdata/ub7+gdata/rt52+gdata/dx2+gdata/lo70+gdata/hh5+gdata/oi10+gdata/fs38

# Load module, always specify version number.
module load R/4.3.1
module load cdo/2.4.3

# Must include `#PBS -l storage=scratch/ab12+gdata/yz98` if the job
# needs access to `/scratch/ab12/` and `/g/data/yz98/`. Details on:
# https://opus.nci.org.au/display/Help/PBS+Directives+Explained

# Run R application
Rscript cmip-damip/compute-sam-cmip6.R