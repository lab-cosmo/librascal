
How to
======

How to edit the documentation
-----------------------------

The format is restructured text. Restructured text work a lot like markdown. However, the syntax is not the same. You can choose whatever editor you prefer to edit the .rst file. Once you edited them save them and remember to create a link to the file you want them to appear.

For example, if you want to add a file describing how_to_do_something to the index, create the how_to_do_someth
ing.rst file and then add 'how_to_do_something' to the index.rst file and it will appear in the HTML.

For a quick reference, here is the cheatsheet for .rst:
https://github.com/ralsina/rst-cheatsheet/blob/master/rst-cheatsheet.rst

How to edit a Readme
--------------------

The format of the README file for git is markdown. Here is a cheat-sheet summarizing how you can use it:
https://guides.github.com/features/mastering-markdown/


How to add a new Machine Learning Method
----------------------------------------

If you want to add a new Machine Learning Method (or in general a method that can be used in Proteus), please, first of all, be aware of the :ref:`code_structure`. Once you know how Proteus is subdivided, you can start thinking about the implementation of the new functionality.

First of all, think about what you need and check if is already implemented (for example, if you need a special C++ structure or if you can use one of the basic types already implemented).

When you know what is necessary for your method, you can proceed in the implementation.

The new method should be written in C++ in the **src** folder and should be named in a way that is clear what the method does (no method_3_bis.cc etc..). The file should start with a header reporting a few information such as the author, date, etc.. You can check the other .cc files in the directory for a quick reference.


Write the C++ method
^^^^^^^^^^^^^^^^^^^^

At this point, you can code the method in the way that you want, but we remind you that it should be included in the ``namespace Proteus {...}``. So for example, your file(s) ``mymethod.h`` will contain:

.. code-block:: c

    #include <wathever_is_needed>
    
    namespace Proteus {
        /**
        * Remember to put comments in a form that Doxygen 
        * so that is clear what they do 
        */
        void function(type& args){
            ...
        }
    }

Once that the method that you need as been coded, you can save the file(s) and proceed. There is no need to change the CMake file. 

Write the binding to the C++ method
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The second task that you need to accomplish is to create a python binding to the C++ subroutine(s) that you just implemented. To accomplish this, you will need to create a ``bind_py_my_method.cc`` file in the **bindings** directory that *at least* the following instruction:
 
.. code-block:: c

    #include <pybind11/pybind11.h>
    #include "mymethod.h"

    using namespace Proteus;
    namespace py=pybind11;
    using namespace pybind11::literals;
    
    /**
    * This create a void function that defines the name of the member
    * that will be called from python, and bind it to function.
    */

    void my_method_binding(py::module& m){
        m.def("name_method",&function);
        /**
        * Add other functions or members that you want to 
        * expose to python.
        * Check the Pybind11 readthedocs for more info
        */
    }

Then, it will be necessary to modify the file ``bind_py_module.cc`` and add your method to those already existing:


.. code-block:: c

    #include <pybind11/pybind11.h>
    #include "mymethod.h"

    using namespace Proteus;
    namespace py=pybind11;
    using namespace pybind11::literals;
    
    /**
    * This is declaring a function of a previous method.
    */
    extern void previous_method(py::module&);

    /**
    * This is declaring the function of your method.
    * This is the line you need to add. 
    */
    extern void my_method(py::module&);

    
    /**
    * This command will expose all the method declared to python,
    * so that it will be possible to import Proteus and use
    * Proteus.previous_method(args)
    * or
    * Proteus.my_method(args)
    */
    PYBIND11_MODULE(_proteus, mod) {
        mod.doc() = "Hello, World!"; //! This is printing the doc.
        previoud_method(mod);
        my_method(mod); //! you also need to add this line
    }


