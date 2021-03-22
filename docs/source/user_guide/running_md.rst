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

The i-PI universal force engine (`ipi-code.org <https://ipi-code.org>`) is a
powerful and flexible simulation engine, capable of many advanced simulation
techniques such as REMD and PIMD and parallel execution of multiple system
replicas.  It is possible to use ``librascal`` as a "force driver" within i-PI's
socket calculator interface: Use `this driver script (on github)
<https://github.com/cosmo-epfl/i-pi/blob/feat/librascal/drivers/py/pes/rascal.py>`
to initialize a saved ``librascal`` model, then run i-PI to connect to the
socket and run dynamics.  Please consult the i-PI documentation for further
information on setting up a simulation and using socket calculators.

ASE (Python)
------------

LAMMPS (under construction)
---------------------------

Others
------

(point to generic MD interface)
