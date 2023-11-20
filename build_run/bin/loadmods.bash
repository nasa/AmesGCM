#!/bin/bash

#source ~/.bashrc
module purge
module load comp-intel
module load mpi-hpe/mpt
module load hdf4/4.2.12
module load hdf5/1.8.18_mpt
module load netcdf/4.4.1.1_mpt
module load szip/2.1.1
#    The following module is needed for nco operators
module load nco/4.6.7