#!/bin/bash

# remove the driver file if it already exists
rm /tmp/ipi_zundel
# make sure rascal can be imported if not installed
export PYTHONPATH="../../../build_b/:$PYTHONPATH"
# path to the i-Pi driver
RASCAL_DRIVER="../../../../i-pi/drivers/py/driver.py"
# i-Pi executable
IPI="i-pi"

# initialize the socket and set up the simulation
$IPI input.xml & sleep 2
# send simulation
$RASCAL_DRIVER -u -a zundel -m rascal -o zundel_model.json,h5o2+.xyz &
wait
