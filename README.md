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
cmake -DCMAKE_BUILD_TYPE=debug  -DENABLE_DOC=ON -DBUILD_TESTS=ON ..
make all
ctest -V
```

cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo  -DBUILD_TESTS=ON ..

* Special flags:
  + INSTALL:
    + empty (default) -> install rascal in the build folder
    + local -> install rascal in the site-package folder of the found python binary
    + pydevelop -> install rascal in the librascal/rascal/lib folder

* Common cmake flag:
  + -DCMAKE_C_COMPILER
  + -DINSTALL
  + -DCMAKE_BUILD_TYPE
  + -DENABLE_DOC
  + -DBUILD_TESTS

RelWithDebInfo

To remove all the cmake files/folders except for the external library (enable glob and remove):
```
shopt -s extglob
rm -fr -- !(external|third-party) 
```
In order to have the readthedocs.org theme for the documentation, please install the following python package:
```Shell
pip install sphinx_rtd_theme
```

To build libRascal with as docker environement:
```
sudo docker build -t test -f ./docker/install_env.dockerfile  .
sudo docker run -it -v /path/to/repo/:/home/user/  test
```
And then follow the instruction in BOOST.md for compilation with boost from conda 


TILL:
Management of derivative relations for fields
"Functional dependency" management to obtain automatically derivatives with chain rule 
Federico:
try to enable CI on github with travis (finalize, enable and think what to do with CI)
Felix:
merge hackaton to master but don't remove the branch
install target (TBD with Till)
start pulling the reference implementation in Python (with Andrea G)
Michele:
Write a whitepaper section
Chiheb:
clean-up and update the tutorial
Everyone:
Make sure new and existing doxygen documentation refers correctly to Order and Layer rather than to Level and Depth
*EFFICIENCY OF BULK KERNEL EVALUATION*
typically we will do operations like diag (AB) where A and B are matrices of the order of 10'000x10'000


TODO:

have NL tests (1-3) that really test all possible features of the implementation