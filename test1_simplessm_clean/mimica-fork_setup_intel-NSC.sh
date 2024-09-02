#!/bin/bash
# Load modules and variables to built mimica

# modules
module use /proj/bolinc/shared/software/modules
module load buildenv-intel/.2018.u3-bare
module load netCDF-HDF5/4.9.2-1.12.2-hpc1
module load FFTW/3.3.8-hpc1
module load CDO/2.3.0-eccodes-aec-cmor-hpc1-intel-2023a-eb

# mimica path 
export MIMICA=/proj/bolinc/users/x_haodi/mimicav5_moding29_28+ssmmod
