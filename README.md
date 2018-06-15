# librascal
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

sudo docker build -t test ./docker/
sudo docker run -it -v /local/git/proteus/:/home/  test
export BOOST_INCLUDEDIR=/usr/local/include/boost/
export BOOST_LIBRARYDIR=/usr/local/lib/
export BOOST_INCLUDEDIR=/usr/include/boost/
export BOOST_LIBRARYDIR=/usr/lib/
export Boost_NO_SYSTEM_PATHS=OFF
export BOOST_NO_SYSTEM_PATHS=OFF

export BOOST_ROOT=/home/boost_1_58_0
export BOOST_LIBRARYDIR=/home/boost_1_58_0/bin.v2/libs/test/build/gcc-5.4.0/release/link-static/threading-multi/
export BOOST_INCLUDEDIR=/home/boost_1_58_0/boost/

export BOOST_ROOT=/local/git/proteus/boost_1_58_0
export  BOOST_LIBRARYDIR=/local/git/proteus/boost_1_58_0/bin.v2/libs/test/build/gcc-5.4.0/release/link-static/threading-multi/
export BOOST_INCLUDEDIR=/local/git/proteus/boost_1_58_0/

export BOOST_ROOT=/local/git/proteus/boost_1_59_0
export BOOST_LIBRARYDIR=/local/git/proteus/boost_1_59_0/bin.v2/libs/test/build/gcc-5.4.0/release/link-static/threading-multi/
export BOOST_INCLUDEDIR=/local/git/proteus/boost_1_59_0/boost/

export BOOST_ROOT=/local/git/proteus/boost_1_67_0


./bootstrap.sh --with-libraries=test



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


