
.. _how_to_contribute:

.. contents::
   :local:


How to contribute
-----------------

If you want to add a new  a method that can be used in Rascal, please, first of all, be aware of the :ref:`code_structure`. Once you know how Rascal is subdivided, you can start thinking about the implementation of the new functionality.

First of all, think about what you need and check if is already implemented (for example, if you need a special C++ structure or if you can use one of the basic types already implemented).

When you know what is necessary for your method, you can proceed in the implementation.

The new method should be written in C++14 in the **src** folder and should be named in a way that is clear what the method does (no method_3_bis.cc etc..). The file should start with a header reporting a few information such as the author, date, etc.. You can check the other .cc files in the directory for a quick reference.
Before, jumping into the next paragraphs, please respect our :ref:`coding convention <coding-convention>`.

Write the C++ method
^^^^^^^^^^^^^^^^^^^^

At this point, you can code the method in the way that you want, but we remind you that it should be included in the ``namespace rascal {...}``. So for example, your file(s) ``mymethod.hh`` will contain:

.. code-block:: c++

  #include <wathever_is_needed>

  namespace rascal {
      /**
      * Remember to put comments in a form that Doxygen
      * so that is clear what they do
      */
      void function(type& args){
          ...
      }
  }

Once that the method that you need as been coded, you can save the file(s) and proceed. Note that adding a ``.cc`` file requires its registration in the folder's CMakeLists.txt so that it can be compiled.
