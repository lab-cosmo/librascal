import json
from itertools import product
from copy import deepcopy
import os
from os.path import join
import sys
import argparse
import numpy as np
from prettyjson import prettyjson

# dump radial and power spectra for methane

root = os.path.abspath("../")
rascal_reference_path = join(root, "reference_data/")
inputs_path = join(rascal_reference_path, "inputs")
read_inputs_path = join("reference_data/", "inputs")
dump_path = join("reference_data/", "tests_only")


def dump_reference_json():
    sys.path.insert(0, join(root, "build/"))
    sys.path.insert(0, join(root, "tests/"))
    from rascal.representations import SphericalInvariants
    from ase.io import read

    np.random.seed(10)
    fns = [
        "diamond_2atom_distorted.json",
        "CaCrP2O7_mvc-11955_symmetrized.json",
        "methane.json",
    ]
    soap_types = ["PowerSpectrum"]
    Nselects = ["all", "all_random", "8_random"]

    sparsification_inputs = []
    for fn, soap_type, Nselect in product(fns, soap_types, Nselects):
        frames = read(join(inputs_path, fn), ":")

        hypers = dict(
            soap_type=soap_type,
            interaction_cutoff=3.5,
            max_radial=2,
            max_angular=2,
            gaussian_sigma_constant=0.4,
            gaussian_sigma_type="Constant",
            cutoff_smooth_width=0.5,
            normalize=False,
            compute_gradients=True,
            expansion_by_species_method="structure wise",
        )

        soap = SphericalInvariants(**hypers)
        managers = soap.transform(frames)
        hyp = deepcopy(hypers)

        # select some features from the possible set
        mapping = soap.get_feature_index_mapping(managers)
        selected_features = {key: [] for key in mapping[0].keys()}
        ids = np.array([key for key in mapping.keys()])
        if Nselect == "all":
            pass
        elif Nselect == "all_random":
            np.random.shuffle(ids)
        elif Nselect == "8_random":
            np.random.shuffle(ids)
            ids = ids[:8]
        else:
            raise NotImplementedError()
        for idx in ids:
            coef_idx = mapping[idx]
            for key in selected_features.keys():
                selected_features[key].append(int(coef_idx[key]))
        # selected_features_global_ids is important for the tests
        selected_features["selected_features_global_ids"] = ids.tolist()
        mapp = dict(coefficient_subselection=selected_features)

        hyp.update(mapp)

        soap_s = SphericalInvariants(**hyp)
        managers_s = soap_s.transform(frames)

        sparsification_inputs.append(
            dict(
                hypers=dict(
                    rep=soap.hypers,
                    rep_sparse=soap_s.hypers,
                    adaptors=json.loads(managers_s.managers.get_parameters()),
                ),
                filename=join(read_inputs_path, fn),
                Nselect=Nselect,
            )
        )

    fn_out = join(root, dump_path, "sparsification_inputs.json")
    print(fn_out)
    with open(fn_out, "w") as f:
        sparsification_inputs_pretty = prettyjson(
            sparsification_inputs, indent=2, maxlinelength=80
        )
        f.write(sparsification_inputs_pretty)


##########################################################################################
##########################################################################################


def main(json_dump):
    if json_dump == True:
        dump_reference_json()


##########################################################################################
##########################################################################################


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-json_dump", action="store_true", help="Switch for dumping json"
    )

    args = parser.parse_args()
    main(args.json_dump)
