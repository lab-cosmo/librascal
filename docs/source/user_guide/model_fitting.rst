.. _model_fitting:

Fitting a Model with Librascal
==============================

Note: This section is under construction.

One of the main things you will want to do with ``librascal``, besides computing
descriptors, is fitting an ML model to a database of structural properties (such
as energies and forces or chemical shifts).  The example notebook
:file:`examples/MLIP_example.ipynb` contains an example of fitting and saving a
model of a potential energy surface to an existing training database.

If you have fit a model to energies and/or forces, then you have a model of the
potential energy surface that you then can use to run molecular dynamics (MD)
simulations.  The notebook :file:`examples/zundel_i-PI.ipynb` has an example of
fitting and running a potential using the included interface to the
`i-PI code (outbound) <http://ipi-code.org>`_.

Finally, you may wish to optimize the speed of the potential by using feature
selection.  This is implemented using the :ref:`filters`
classes and demonstrated in the notebook :file:`examples/Feature_selection_example.ipynb`.
