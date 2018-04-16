# proteus
A scalable and versatile fingerprint and machine-learning code

How to install:
* Need to have the programs git, gcc (or other c++14 compiler), boost (unit_test_framework), doxygen
* Need the python packages Sphinx, Breathe 
* To build the program: 
```Shell
mkdir build 
cd build 
cmake .. 
make
``` 
* To make development documentation: first enable the documentation building with ccmake, then
```Shell
cd build 
cmake ..
make dev_doc
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
Chiheb:
Start writing down the blueprints from the google doc to the docs
*EFFICIENCY OF BULK KERNEL EVALUATION*
typically we will do operations like diag (AB) where A and B are matrices of the order of 10'000x10'000


