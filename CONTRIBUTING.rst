Contributing to librascal
********************************************************************************

Thank you for contributing to libRascal! To make sure your pull request gets
accepted with a minimum of hassle for everyone involved, please make sure first
and foremost that you've read the ``README`` (You'd think this would be obvious,
but far too many people ignore ``README``\ s in practice). It is also a good
idea to browse the `developer's guide
<https://cosmo-epfl.github.io/librascal/dev_guide/developer.html>`_, especially
the `coding conventions
<https://cosmo-epfl.github.io/librascal/dev_guide/coding-convention.html>`_ for
C++ code (Python code generally follows PEP 8).


Development dependencies
================================================================================

All python packages required for development (including all python dependencies
of the librascal project) can be installed with

.. code:: shell

   pip install -r requirements/dev.txt


Summary of the review process
================================================================================

We have a pull request template; please use it (especially for new features)! In
particular, make sure you've checked off all the items listed under "Tasks
before review". If you haven't been able to complete all these tasks, then your
code isn't ready for review; consider submitting a "Draft pull request" instead
(or waiting to request reviewers).

A detailed description of the review process can be found in a separate document
about the `review process
<https://cosmo-epfl.github.io/librascal/dev_guide/review_process.html>`_. Here
is a short summary of the important points:

For developers

 * We want a clean and tested master branch, which is why we review code and use
   continuous integration tools.
 * Contact us if you have any questions.
 * Start the discussion about your contribution early by using the
   Draft-Pull-Request feature of Github.
 * Use our Pull-Request template in preparing your PR.
 * Provide tests for the new/changed functionalities in the test suite of the
   library.
 * Make sure your proposed changes pass the all existing tests, change (not
   deactivate) them if necessary.
 * Document your contribution.

For reviewers

 * The first person to make a full review is responsible for seeing it through.
 * Single line comments are ok, but consider a proper review.
 * If you started a review, you are responsible for seeing it through.
 * Finish the PR fast by reviewing the requested changes.
 * If you start changing code on a PR as a reviewer, ask another person to
   review you code.

Jupyter Notebooks
================================================================================


If you do have to commit notebooks to version control, make sure you've
installed the ``nbstripout`` extension (available on `github
<https://github.com/kynan/nbstripout#installation>`_ and `PyPI
<https://pypi.org/project/nbstripout/>`_) and that it works properly (i.e.
you're not committing the output of notebook cells). Any pull requests that
include diffs of notebook outputs will be sent back without further comments.
The two exceptions to this rule are:

1. Deleting existing, committed outputs (a rebase may be requested to prevent
   those outputs from making it into the history in the first place)

2. Stable outputs for example notebooks meant to be processed to HTML for the
   `tutorials page
   <https://cosmo-epfl.github.io/librascal/tutorials/tutorials.html>`_. "Stable"
   here means these are the final outputs meant to be shown to the public and
   won't be changed unless errors or omissions are discovered, or the tutorial
   is later updated or expanded. Note that tutorials that rely on volatile code
   features should *not* have their outputs included, since in that case the
   tutorials themselves aren't stable; instead, you should use the
   auto-execution feature of ``nbsphinx`` (see `here
   <https://nbsphinx.readthedocs.io/en/latest/executing-notebooks.html>`_ for
   details).

Nonetheless, it is highly discouraged to contribute code in the form of
notebooks; even with filters like ``nbstripout`` they're a hassle to use in
version control. Use them only for comprehensive tutorials or *stable* examples
that are either meant to be run *interactively* or are meant to be processed by
`sphinx` (`nbsphinx <https://nbsphinx.readthedocs.io/en/latest/>`_) for
inclusion in the `introductive examples
<https://cosmo-epfl.github.io/librascal/examples/examples.html>`_.

After installing ``nbstripout``, activate it for this project by running:

.. code:: shell

   nbstripout --install --attributes .gitattributes

from the top-level repository directory. Please note that that ``nbstripout``
will not strip output from cells with the metadata fields ``keep_output`` or
``init_cell`` set to ``True``, so use these fields judiciously. You can ignore
these settings with the following command:

.. code:: shell

   git config filter.nbstripout.extrakeys '\
      cell.metadata.keep_output cell.metadata.init_cell'

(The keys ``metadata.kernel_spec.name`` and
``metadata.kernel_spec.display_name`` may also be useful to reduce diff noise.)

Note also that ``nbstripout`` will not strip output from cells with the metadata
fields ``keep_output`` or ``init_cell`` set to ``True``. Please use these fields
judiciously.


Developer tools
================================================================================

Tools to make life simpler for developers

Deepclean
--------------------------------------------------------------------------------

To remove all the cmake files/folders except for the external libraries:

.. code:: shell

   make deepclean

Automatic code formatting
--------------------------------------------------------------------------------

To help developers conform their contribution to the coding convention, the
formatting of new functionalities can be automated using clang-format (for the
c++ files) and black (for the python files). The .clang-format and .pycodestyle
files define common settings to be used.

To enable these functionalities (optional) you can install these tools with:

.. code:: shell

   sudo apt-get install clang-format-8
   pip install black

The automatic formatting of the c++ and python files can be triggered by:

.. code:: shell

   cd build
   cmake ..
   make pretty-cpp
   make pretty-python


Note that the use of auto-formatters isn't mandatory, but your code *must* pass
the various linters that are run with ``make lint``; this is part of the
automated CI build. Any code that has been auto-formatted should pass the
linters; please open an issue if this is not the case for your code. If the
linter fails on your code, the pull request will be sent back for revision
without any further comments.

Please use these tools with caution as they can potentially introduce unwanted
changes to the code. If code needs to be specifically excluded from auto
formatting, e.g. a matrix which should be human-readable, code comments tells
the formatters to ignore lines:

- C++

  .. code:: C++

     // clang-format off
     SOME CODE TO IGNORE
     // clang-format on

- python

  .. code:: python

     SOME LINE TO IGNORE # noqa

  where ``noqa`` stands for ``no`` ``q``\ uality ``a``\ ssurance.

