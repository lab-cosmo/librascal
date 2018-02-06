cmake_minimum_required(VERSION 2.8.12)
project(proteus)

add_subdirectory(pybind11)
pybind11_add_module(proteus src/module.cc)
#add_subdirectory(proteus)