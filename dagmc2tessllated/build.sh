#!/bin/bash

g++ test.cc -I$HOME/dev/dagmc/include -I$HOME/opt/moab/include -I$HOME/opt/geant4/include/Geant4/ -L$HOME/dev/dagmc/lib -ldagsolid -ldagmc -L$HOME/opt/moab/lib -lMOAB -L$HOME/opt/geant4/lib/ -lG4geometry -lG4clhep -o test

g++ test_one.cc -I$HOME/dev/dagmc/include -I$HOME/opt/moab/include -I$HOME/opt/geant4/include/Geant4/ -L$HOME/dev/dagmc/lib -ldagsolid -ldagmc -L$HOME/opt/moab/lib -lMOAB -L$HOME/opt/geant4/lib/ -lG4geometry -lG4clhep -o test_one
