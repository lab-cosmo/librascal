.. _running_md:

Running MD Simulations with a Librascal Model
=============================================

Once you have fitted a model with ``librascal`` (see :ref:`model_fitting`), you will
probably want to use it to do something interesting -- for example, running a
molecular dynamics (MD) simulation.  Librascal is not an MD code, so it provides
interfaces to other codes that specialize in running molecular simulations;
``librascal`` only calculates descriptors (and usually also energies and forces)
given an atomic configuration.  Below is a list of codes that can currently use
``librascal`` to run MD simulations.

i-PI
----

The i-PI universal force engine (`ipi-code.org <https://ipi-code.org>`_) is a
powerful and flexible simulation engine, capable of many advanced simulation
techniques such as REMD and PIMD and parallel execution of multiple system
replicas.  It is possible to use ``librascal`` as a "force driver" within i-PI's
socket calculator interface: Use `this driver script (on github)
<https://github.com/cosmo-epfl/i-pi/blob/feat/librascal/drivers/py/pes/rascal.py>`_
to initialize a saved ``librascal`` model, then run i-PI to connect to the
socket and run dynamics.  Please consult the i-PI documentation for further
information on setting up a simulation and using socket calculators.

ASE (Python)
------------

Librascal also provides a Python Calculator compatible with the `ASE
<https://wiki.fysik.dtu.dk/ase/>`_ package, which is capable of running serial
MD (``librascal`` does not yet parallelize the energy or force calculation).
The Calculator can be used for parallel evaluation of energies/forces for a set
of *existing* structures with the help of a standard Python parallelism tool,
such as `ipyparallel <https://ipyparallel.readthedocs.io/en/latest/>`_.

LAMMPS (under construction)
---------------------------

`LAMMPS <https://lammps.sandia.gov/>`_ is one of the most widely used MD codes
and supports a wide array of energy and force models (including several machine
learning potentials), as well as parallel MD via spatial domain decomposition.
An interface to LAMMPS is currently in the works; stay tuned for more!

Others
------

If you would like ``librascal`` to support another MD code, please let us know!
The best way to do this is `open an issue on our github
<https://github.com/cosmo-epfl/librascal/issues>`_ with the ``enhancement``
label.

In the meantime, if you can make your MD code call Python code, you may be able
to use the :class:`.GenericMDCalculator` class defined in
`bindings/rascal/models/IP_generic_md.py` -- which can be used as-is or as the
starting point for a more specialized, optimized interface.
