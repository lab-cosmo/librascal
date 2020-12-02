import collections
import os
import sys

import ase.io
import numpy as np

import rascal.utils
import rascal.utils.io
import rascal.models
from rascal.models import gaptools
from tqdm import tqdm


WORKDIR = "potential_default"


# TODO do_gradients is assumed; make a version without forces
#      (and without energies)
# TODO parameter parsing and validation is sorta entwined with the work here;
#      consider offloading it to something like argparse
def fit_save_model(parameters):
    WORKDIR = parameters["working_directory"]
    geom_subset = parameters.get("structure_subset", ":")
    source_geoms = ase.io.read(parameters["structure_filename"], geom_subset)
    rep, soaps, sparse_points = gaptools.calculate_and_sparsify(
        tqdm(source_geoms), parameters["soap_hypers"], parameters["n_sparse"]
    )
    (kobj, kernel_sparse, kernel_energies, kernel_forces) = gaptools.compute_kernels(
        rep, soaps, sparse_points, parameters.get("soap_power", 2)
    )
    eparam_name = parameters.get("energy_parameter_name", "energy")
    fparam_name = parameters.get("force_parameter_name", "forces")
    energies = np.array([geom.info[eparam_name] for geom in source_geoms])
    forces = np.array([geom.arrays[fparam_name] for geom in source_geoms])
    energy_baseline_in = parameters["atom_energy_baseline"]
    if isinstance(energy_baseline_in, collections.abc.Mapping):
        energy_baseline = energy_baseline_in
    else:
        all_species = set()
        for geom in source_geoms:
            all_species = all_species.union(geom.get_atomic_numbers())
            # Convert to int because otherwise it's a numpy type
            # (which doesn't play well with json)
            all_species = set(int(sp) for sp in all_species)
        energy_baseline = {species: energy_baseline_in for species in all_species}
    energy_delta = parameters.get("energy_delta", 1.0)
    energy_regularizer = parameters["energy_regularizer"] / energy_delta
    force_regularizer = parameters["force_regularizer"] / energy_delta
    weights = gaptools.fit_gap_simple(
        source_geoms,
        kernel_sparse,
        energies,
        kernel_energies,
        energy_regularizer,
        energy_baseline,
        forces,
        kernel_forces,
        force_regularizer,
    )
    np.save(os.path.join(WORKDIR, "weights"), weights)
    model = rascal.models.KRR(weights, kobj, sparse_points, energy_baseline)
    rascal.utils.dump_obj(os.path.join(WORKDIR, "gap_model.json"), model)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise RuntimeError(
            "Must provide the name of a JSON parameter file"
            "as the first argument (see the examples/ directory"
            " for, well, an example)"
        )
    elif len(sys.argv) > 2:
        raise RuntimeError("Too many arguments")
    else:
        fit_save_model(rascal.utils.io.load_json(sys.argv[1]))
