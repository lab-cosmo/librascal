Contributing to libRascal
-------------------------

Thank you for contributing to libRascal!  To make sure your contribution gets
accepted with a minimum of hassle for everyone involved, please make sure first
and foremost that you've read the `README`, especially the section entitled
"Helpers for Developers", which contains details on the formatting functions
that you should run on your code before sending a pull request for review.
(Note that the use of auto-formatters isn't mandatory, but your code *must* pass
the various linters that are run with `make lint`; this is part of the automated
CI build.  Any code that has been auto-formatted should pass the linters; please
open an issue if this is not the case for your code).  If the linter fails on
your code, the pull request will be sent back for revision without any further
comments.

It is also a good idea to browse the `developer's guide
<https://cosmo-epfl.github.io/librascal/dev_guide/developer.html>`_, especially
the `coding conventions
<https://cosmo-epfl.github.io/librascal/dev_guide/coding-convention.html>`_ for
C++ code (Python code generally follows PEP 8).

Jupyter Notebooks
=================

Please read the README section on this topic.  In general, avoid contributing
code in notebooks wherever possible.  If you do have to commit notebooks to
version control, make sure you've installed the `nbstripout` extension and that
it works properly (i.e. you're not committing the output of notebook cells).
