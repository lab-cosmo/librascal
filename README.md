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
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

* Python requirements: python3.5 and newer, numpy, scipy, ASE (https://wiki.fysik.dtu.dk/ase/index.html), cpplint(optional), Sphinx(optional). To install these packages you could run:
```Shell
pip install numpy scipy ase cpplint sphinx sphinx_rtd_theme
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
cmake -DCMAKE_BUILD_TYPE=Debug  -DENABLE_DOC=ON -DBUILD_TESTS=ON ..
make all
ctest -V
```

* To check for conformity with the c++ code convention:
```Shell
cd build
cmake -DCMAKE_BUILD_TYPE=Debug  -DENABLE_DOC=ON -DBUILD_TESTS=ON ..
make lint
```

* To build with optimization and debug info:
```Shell
cd build
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo -DBUILD_TESTS=ON  CMAKE_C_FLAGS_RELWITHDEBUBINFO="-03 -g -DNDEBUG" ..
make -j 4
ctest -V
```

* To help developers conform their contribution to the coding convention, the formating of new functionalities can be automated using clang-format (for the c++ files) and autopep8 (for the python files). The .clang-format and .pycodestyle files define common settings to be used. To enable these functionalities (optional) you can install these tools with:
```Shell
sudo apt-get install clang-format
pip install autopep8
```
* The automatic formating of the c++ and python files can be trigered by:
```Shell
cd build
cmake ..
make pretty-cpp
make pretty-python
```
Please use these tools with caution as they can potentially introduce unwanted changes to the code.
If code needs to be specifically excluded from auto formatting, e.g. a matrix which should be human-readable, code comments tells the formatters to ignore lines:

C++
```
// clang-format off
SOME CODE TO IGNORE
// clang-format on
```

python
```
SOME LINE TO IGNORE # noqa
```
where `noqa` stands for `no` `q`uality `a`ssurance.

* Common cmake flag:
  + -DCMAKE_C_COMPILER
  + -DBUILD_BINDINGS
  + -DUSER
  + -DINSTALL_PATH
  + -DCMAKE_BUILD_TYPE
  + -DENABLE_DOC
  + -DBUILD_TESTS

* Special flags:
  + -DBUILD_BINDINGS:
    + ON (default) -> build python binding
    + OFF -> does not build python binding
  + -DINSTALL_PATH:
    + empty (default) -> does not install in a custom folder
    + custom string -> root path for the installation
  + -DUSER:
    + OFF (default) -> changes nothing
    + ON -> install root is in the user's home directory, i.e. ~/.local/

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
