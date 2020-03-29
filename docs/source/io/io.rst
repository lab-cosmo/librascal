.. _io:

Saving and Loading in Librascal
===============================

The library aims at being able to load/save three types of objects: atomic structures,
their computed representation and trained machine learning models.
Many tools already exist for reading atomic srtuctures in python and c++ so we rely on
ASE_ on the python side and
support the ASE `.json` atomic structure file format on the c++ side
(:cpp:class:`AtomicStructure <rascal::AtomicStructure>`). An example is available
in `json_structure.cc` and `spherical_invariants_example.cc`.
At the moment we do not support saving representations in native Librascal format
but on the python side they can be converted to numpy arrays with the `get_features()`
function and then saved.
The serialization of Librascal's native c++ objects rely on the `json type
<https://github.com/nlohmann/json>`_
but only a few object have this capability (at the moment only the ones required for saving a trained KRR model).
On the other hand Librascal's python objects that inherit from the :py:class:`rascal.utils.BaseIO` can be serialized
to a dictionary (see :py:meth:`rascal.utils.BaseIO.to_dict` function) or a file, :py:meth:`rascal.utils.BaseIO.dump_obj`
which supports the `JSON` format. Note that numpy arrays larger than 200MB will be saved in a separate `.npy` file
referenced in the main object file (both files should be in the same folder).
The function :py:meth:`rascal.utils.BaseIO.load_obj` loads a Librascal python object.
An example is available in the `MLIP_example.ipynb` jupyter notebook.

.. _ASE: https://wiki.fysik.dtu.dk/ase/ase/io/io.html
