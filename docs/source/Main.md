	\section A The main idea
This is the documentation of Rascal <a href="https://github.com/cosmo-epfl/librascal"><b>(GitHub)</b></a> - free scalable and versatile library to generate representations for atomic-scale learning. It collects several different algorithms that can be used to create fingerprints, perform dimensionality reduction and fit, atomistic and finite element calculations. Rascal calculations is based on the many papers, describing such things as molecular descriptors and the Gaussian process regression in application to molecular properties calculations. Rascal was developed by the team of 
<a href="https://cosmo.epfl.ch/">Laboratory of Computational Science and Modelling</a> in EPFL. 

Rascal as of now is thought as standalone code. However, we aim to provide enough flexibility to interface it with other codes such as LAMMPS and PLUMED-2.0. It can be used as a C++ library as well as a python module. To be able to call it from python, we have used the pybind11 library.

Although at the moment is a serial-only code, we aim to write it in MPI so that it will be possible to take advantage of parallelization to speed up the calculations significantly.

It comes with a GNU Lesser General Public License of version 3, which means that it can be modified and freely distribute, although we take no responsibility for its misuse.

	\section B Code structure
The code is divided mainly in two parts: a pure C++ part and a python-binding interface.

The subroutines of the code that performs the expensive part of the calculation are written in C++14, and are collected in the **/src/** directory, as in the case of the :cpp:func:`cdist <cdist>`. This folder is completely agnostic of the python binding, and it should be kept in this way.

The python-bindings is obtained through Pybind11, and the binding subroutines are included in the **/bindings/** folder. Here, the bind_py_module.cc is the file that contains the main binding of the C++ package (in other words, is what is needed to use the syntax ``import Rascal``). The binding for each submodule of Rascal and its members are collected in files named bind_py_METHOD.cc. We decided to employ Pybind11 because of its seamless integration between Eigen and numpy. This allows the developer (and the user) to code fast and efficient algorithms easily, without losing the power of the C++ linear algebra as well as numpy simplicity. For more reference, please consult `Eigen interface from pybind11 <http://pybind11.readthedocs.io/en/stable/advanced/cast/eigen.html?highlight=eigen#pass-by-reference>`_.

	\section C Getting started
	To see how to install Rascal, follow the <a href="md_Installation.html"><b> installation link</b></a>.

	\section D Tutorial
	As a tutorial, you can have a look at two python notebooks to see how you can proceed from the dataset to molecular properties. In case you are interested in that, follow the <a href="usergroup0.html"><b>tutorial link</b></a>.

	\section E Theory
	To learn more about the theory standing behind the Rascal, follow the <a href="SOAP_documentation.pdf"><b>review link</b></a>. It contains the paper describing the math used in calculations with some references. For [kind of] full reference list click on <a href="md_References.html"><b>references link</b></a>.
	\section F Developers guide
	If you want to contribute to project, you can follow the <a href="usergroup2.html"><b>guide link</b></a> and see the tips for that.
	\section G Auto documentation
	To learn more about the code structure, you can take a look at the auto-generated documentation of the code: <br>
<a href="namespaces.html"><b>Namespaces</b></a> <br>
<a href="annotated.html"><b>Classes</b></a> <br>
<a href="files.html"><b>Files</b></a> 
