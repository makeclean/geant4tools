#!/bin/bash

MOAB_ROOT=/home/andrewdavisdavis/opt/moab
DAGMC_ROOT=/home/andrewdavisdavis/dev/dagmc
HDF5_ROOT=/usr/lib/x86_64-linux-gnu/hdf5/serial

g++ -g -std=c++11 main.cc -I$MOAB_ROOT/include -I$DAGMC_ROOT/include -I$HDF5_ROOT/include -L$MOAB_ROOT/lib -lMOAB -L$DAGMC_ROOT/lib -ldagmc -lpyne_dagmc -L$HDF5_ROOT/lib -lhdf5_serial -lxerces-c -o main
