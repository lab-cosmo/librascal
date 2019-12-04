# Rascal Scripts

This folder contains some utility scripts used by active developers of the library.
They are meant for development purposes only.

## generate_*.py

These scripts use the library and its python bindings to produce reference outputs used in the c++ regression tests.
To use them, compile the library with the bindings and make sure that the state of the code will result in an output corresponding to a valid reference !

## developer_utils.py

Used by the `make` taget `deepclean` to remove all files from the current build folder except for the downloaded libraries.

## fixheader.vim and license-header.cpp.txt

`vim` script to add a licence to all c++ sources of the library.
