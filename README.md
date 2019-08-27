# librascal
A scalable and versatile fingerprint and machine-learning code

## How to install:
* Need to have the programs git, gcc/clang (with c++14 support), cmake and python 3.6 with numpy, scipy and ASE (https://wiki.fysik.dtu.dk/ase/index.html). On Ubuntu run:
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


## Alternative build:

The library supports several alternative build that have additional dependencies.
Note that the `ncurses` GUI for cmake (ccmake) is quite helpful to customize the build options.

1. Tests

    Librascal source code is extensively tested (both c++ and python). The BOOST unit_test_framework is requiered to build the tests (see BOOST.md for further details on how to install the boost library).
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

2. Build Type

    Several build types are available Release (default), Debug and RelWithDebInfo. To build an alternative mode
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

3. Documentation

    The documentation relies on the sphinx, breathe, doxygen and graphviz packages. To install them on ubuntu:
    ```Shell
    pip3 install sphinx sphinx_rtd_theme breathe
    sudo apt-get install doxygen graphviz
    ```
    Then to build the documentation run:
    ```Shell
    cd build
    cmake -DENABLE_DOC=ON  ..
    make dev_doc
    ```

4. Helpers for Developers

     * To remove all the cmake files/folders except for the external library (enable glob and remove):
    ```Shell
    shopt -s extglob
    rm -fr -- !(external|third-party)
    ```
    * To help developers conform their contribution to the coding convention, the formatting of new functionalities can be automated using clang-format (for the c++ files) and autopep8 (for the python files). The .clang-format and .pycodestyle files define common settings to be used.

      To enable these functionalities (optional) you can install these tools with:
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

5. Bindings

    Librascal relies on the pybind11 library to automate the generation of the python bindings which are built by default. Nevertheless, to build only the c++ library:
    ```Shell
    cd build
    cmake -DBUILD_BINDINGS=OFF ..
    make
    ```

## Miscellaneous Information

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


To build librascal as a docker environment:
```
sudo docker build -t test -f ./docker/install_env.dockerfile  .
sudo docker run -it -v /path/to/repo/:/home/user/  test
```
