# librascal
A scalable and versatile fingerprint and machine-learning code

How to install:
* Need to have the programs git, gcc (or other c++14 compiler), boost (unit_test_framework, see BOOST.md for further details on how to install the boost library), doxygen
* Need the python packages Sphinx, Breathe 
* You can use either cmake or (recommended) ccmake to configure the Makefile. If using ccmake, you should press *c* one or more times until the option to generate the configuration appears, then type *g*.
* To build the program: 
```Shell
mkdir build 
cd build 
ccmake .. 
make
``` 
* To make development documentation: first enable the documentation building with ccmake, then
```Shell
cd build 
ccmake ..
make dev_doc
``` 

* To build for development:
```Shell
cd build 
cmake -DCMAKE_BUILD_TYPE=Debug  -DENABLE_DOC=ON -DBUILD_TESTS=ON -DCMAKE_C_COMPILER=/usr/local/bin/gcc-7 -DCMAKE_CXX_COMPILER=/usr/local/bin/g++-7 ..
make all
ctest -V
```
To remove all the cmake files/folders except for the external library (enable glob and remove):
```
shopt -s extglob
rm -fr -- !(external|third-party) 
```
In order to have the readthedocs.org theme for the documentation, please install the following python package:
```Shell
pip install sphinx_rtd_theme
```


TILL:
Test + compilation of documentation in cmake + how to add a test
Make sure that one can compile on fidis / daint without docs and tests and it will work based on existing libraries
Libboost-test dependency 
Stub for the neighbor list
Management of derivative relations for fields
"Functional dependency" management to obtain automatically derivatives with chain rule 
s/Field/Property/g
Federico:
comments to c++ + how to bindings/function/class
try to enable CI on github with travis
Felix:
the rest of the tutorial
start pulling the reference implementation in Python (with Andrea G)
Michele:
Write a whitepaper section
Everyone:
Make sure new and existing doxygen documentation refers correctly to Order and Layer rather than to Level and Depth
*EFFICIENCY OF BULK KERNEL EVALUATION*
typically we will do operations like diag (AB) where A and B are matrices of the order of 10'000x10'000


