How to edit the documentation
=============================

The format is restructured text. Restructured text work a lot like markdown. However, the syntax is not the same. You can choose whatever editor you prefer to edit the .rst file. Once you edited them save them and remember to create a link to the file you want them to appear.

For example, if you want to add a file describing how_to_do_something to the index, create the how_to_do_someth
ing.rst file and then add 'how_to_do_something' to the index.rst file and it will appear in the HTML.

For a quick reference, here is the cheatsheet for .rst:
https://github.com/ralsina/rst-cheatsheet/blob/master/rst-cheatsheet.rst

How to edit a Readme
====================

The format of the README file for git is markdown. Here is a cheat-sheet summarizing how you can use it:
https://guides.github.com/features/mastering-markdown/


How to add a new Machine Learning Method.
=========================================

f you want to add a new Machine Learning Method (or in general a method that can be used in Proteus), please, first of all, be aware of the :ref:`code_structure`. Once you know how Proteus is subdivided, you can start thinking about the implementation of the new functionality.

First of all, think about what you need and check if is already implemented (for example, if you need a special C++ structure or if you can use one of the basic types already implemented). 

When you know what is necessary for your method, you can proceed in the implementation.

The new method should be written in C++ in the **\src\** folder, and should be named in a way that is apparent what the method does (no method_3_bis.cc etc..). The file should starts with a header reporting a few information such as the author, date etc.. You can check the other .cc files in the directory for a quick reference. 

At this point, you can code the method in the way that you want, but we remind you that it should be included in the ``namespace Proteus {...}``.

Once that the method that you need as been coded in C++, you can save you .cc and .h files and proceed. There is no need to change the CMake file.


