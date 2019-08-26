
# Installation
* Need to have the programs Git, gcc/clang (with c++14 support), CMake and Python 3.6 with numpy, scipy and ASE (https://wiki.fysik.dtu.dk/ase/index.html). On Ubuntu run:
```Shell
sudo apt-get install git gcc cmake
pip3 install numpy scipy ase
```
* Then you can create a build folder, generate the make files and build the library:
```Shell
mkdir build
cd build
cmake ..
make
```


## Alternative build

The library supports several alternative build that have additional depencies.
Note that the curse gui for cmake (ccmake) is quite helpfull to customize the build options.

1. <b>Tests</b>

    Librascal source code is extensively tested (both c++ and python). The BOOST unit_test_framework is requiered to build the tests (see below for further details on how to install the boost library).
    To build and run the tests:
```Shell
cd build
cmake -DBUILD_TESTS=ON ..
make
ctest -V
```
    In addition to testing the behaviour of the code, the test suite also check for formatting compliance with the clang-format and autopep8 packages (these dependencies are optional).
    To install these dependencies on ubuntu:
```Shell
sudo apt-get install clang-format
pip3 install autopep8
```

2. <b>Build Type</b>

    Several build types are available Release (default), Debug and RelWithDebInfo. To build an alternative mode:
```Shell
cd build
cmake -DCMAKE_BUILD_TYPE=Debug  ..
make
```
    Or
```Shell
cd build
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo  CMAKE_C_FLAGS_RELWITHDEBUBINFO="-03 -g -DNDEBUG" ..
make
```

3. <b>Documentation</b>

    The current documentation relies on the doxygen package. To install it on ubuntu:
```Shell
sudo apt-get install doxygen
```
    The documentation is located in the librascal/docs/documentation/html folder. The source files for the documentation are located in the librascal/docs/src folder. 
    To rebuild the documentation, run the 
```Shell
doxygen Config
```
    in the librascal/docs/src folder.
4. <b>Bindings</b>

      Librascal relies on the pybind11 library to automate the generation of the python bindings which are built by default. Nevertheless, to build only the c++ library:
```
cd build
cmake -DBUILD_BINDINGS=OFF ..
make
```

5. <b>Helpers for Developers</b>

* To remove all the cmake files/folders except for the external library (enable glob and remove):
```Shell
shopt -s extglob
rm -fr -- !(external|third-party)
```
* To help developers conform their contribution to the coding convention, the formating of new functionalities can be automated using clang-format (for the c++ files) and autopep8 (for the python files). The .clang-format and .pycodestyle files define common settings to be used.

* To enable these functionalities (optional) you can install these tools with:
```Shell
sudo apt-get install clang-format
pip install autopep8
```
      The automatic formating of the c++ and python files can be trigered by:
```Shell
cd build
cmake ..
make pretty-cpp
make pretty-python
```
      Please use these tools with caution as they can potentially introduce unwanted changes to the code.
      If code needs to be specifically excluded from auto formatting, e.g. a matrix which should be human-readable, code comments tells the formatters to ignore lines:

* C++
```Shell
// clang-format off
SOME CODE TO IGNORE
// clang-format on
```

* Python

```Shell
SOME LINE TO IGNORE # noqa
```
      where <b>`noqa`</b> stands for <b>no</b> <b>q</b>uality <b>a</b>ssurance.

## Misceleaneous Information

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


To build libRascal with as docker environement:
```
sudo docker build -t test -f ./docker/install_env.dockerfile  .
sudo docker run -it -v /path/to/repo/:/home/user/  test
```

## How to install Boost unit test framework library

We are using the boost library for our c++ testing framework following the 1.62 and onward syntax. There are several solutions to install it.

## Ubuntu/Debian system wide

```
sudo apt-get install libboost-test-dev
```
Depending on your version of these distributions, you might not get a recent enough verison of the boost library.

## Mac
It seems that the brew recepy of boost is broken (compilation errors). Custom or Conda installations are advised. 

## Custom installation
We use 7z to unpack (sudo apt-get install p7zip-full to install on ubuntu)
```
export BOOST_DIR=/path/to/download/directory
wget https://sourceforge.net/projects/boost/files/boost/1.67.0/boost_1_67_0.7z/download -O $BOOST_DIR/boost_1_67_0.7z
7z x $BOOST_DIR/boost_1_67_0.7z
cd $BOOST_DIR/boost_1_67_0
mkdir boost_root
./bootstrap.sh --with-libraries=test --with-toolset=gcc --prefix=$BOOST_DIR/boost_root 
./b2 install
export Boost_NO_SYSTEM_PATHS=ON
export BOOST_ROOT=$BOOST_DIR/boost_root
```
Cmake FindBoost might not find this installation of Boost Test so we set Boost_NO_SYSTEM_PATHS to avoid finding another installation of Boost and we set BOOST_ROOT to this new installation.

## Using Conda
Conda distributes Boost library in the conda forge channel.
+ First install miniconda if not already available:
```
export MINICONDA_DOWNLOAD_DIR=/path/to/download/directory
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O $MINICONDA_DOWNLOAD_DIR/miniconda.sh
bash $MINICONDA_DOWNLOAD_DIR/miniconda.sh
bash $MINICONDA_DOWNLOAD_DIR/miniconda.sh -b -p $HOME/miniconda
export PATH="$HOME/miniconda/bin:$PATH" 
echo 'export PATH="$HOME/miniconda/bin:$PATH"' >> $HOME/.bashrc
source $HOME/.bashrc
conda create -n myenv python=3.6
```

+ Install boost from conda forge
```
source activate myenv
conda install -c conda-forge boost
export Boost_NO_SYSTEM_PATHS=ON 
export BOOST_ROOT=/path/to/myenv
```
Following the previous example, /path/to/myenv -> $HOME/miniconda/envs/myenv.

Cmake FindBoost might not find this installation of Boost so we set Boost_NO_SYSTEM_PATHS to avoid finding another installation of Boost and we set BOOST_ROOT to this new installation.
```
