#!/bin/bash

# remove the driver file if it already exists
rm /tmp/ipi_zundel
# make sure rascal can be imported if it is not installed
#export PYTHONPATH="../../../build/:$PYTHONPATH"
# source env.sh in the i-pi folder so that the driver is in the path
# source $IPI_PATH/env.sh
RASCAL_DRIVER="i-pi-py_driver"
# i-Pi executable
IPI="i-pi"

# initialize the socket and set up the simulation
$IPI input.xml & sleep 2
# send simulation
$RASCAL_DRIVER -u -a zundel -m rascal -o zundel_model.json,h5o2+.extxyz &
wait
