import unittest
import numpy as np
import ase.io

from rascal.utils import filter, FPSFilter, CURFilter
from rascal.representations import SphericalInvariants as SOAP


class IndexConversionTest(unittest.TestCase):

    """Test index conversion utilities for filtering"""

    def setUp(self):
        self.example_frames = ase.io.read(
            "reference_data/inputs/small_molecules-20.json", ":5"
        )
        global_species = set()
        for frame in self.example_frames:
            global_species.update([int(sp) for sp in frame.get_atomic_numbers()])
        self.global_species = global_species
        self.repr = SOAP(
            max_radial=3,
            max_angular=0,
            soap_type="RadialSpectrum",
            interaction_cutoff=4.0,
            cutoff_smooth_width=1.0,
            gaussian_sigma_type="Constant",
            gaussian_sigma_constant=0.3,
            expansion_by_species_method="user defined",
            global_species=sorted(list(global_species)),
        )
        self.managers = self.repr.transform(self.example_frames)
        self.example_features = self.managers.get_features(self.repr)

    def test_split_by_sp(self):
        """Test the feature matrix species splitter"""
        # Let's try a basic example matrix first
        X = np.concatenate(
            [frame.get_atomic_numbers() for frame in self.example_frames]
        )[:, np.newaxis]
        sps = list(self.global_species)
        X_by_sp = filter._split_feature_matrix_by_species(self.managers, X, sps)
        for sp in sps:
            self.assertTrue(np.all(X_by_sp[sp] == sp))
        # Now with the actual feature matrix
        X_by_sp = filter._split_feature_matrix_by_species(
            self.managers, self.example_features, sps
        )
        X_by_sp_manual = {}
        atoms_species = X.flatten()
        for sp in sps:
            X_by_sp_manual = self.example_features[atoms_species == sp]
            self.assertTrue(np.all(X_by_sp_manual == X_by_sp[sp]))

    def test_indices_global(self):
        """Test the global-to-perstructure index transformation"""
        example_idces_global = [0, 2, 16, 17, 82]
        example_idces_perstructure = filter._indices_manager_to_perstructure(
            self.managers, example_idces_global
        )
        self.assertEqual(example_idces_perstructure, [[0, 2, 16], [0], [], [], [19]])
        # It should also work with an ASE list of atoms
        example_idces_perstructure = filter._indices_manager_to_perstructure(
            self.example_frames, example_idces_global
        )
        self.assertEqual(example_idces_perstructure, [[0, 2, 16], [0], [], [], [19]])

    def test_indices_global_out_of_range(self):
        """Test the index transformer with out-of-range indices"""
        with self.assertRaisesRegex(
            ValueError, rf"Selected index\(es\): \[83\] out of range"
        ):
            filter._indices_manager_to_perstructure(
                self.managers,
                [
                    83,
                ],
            )
        bad_indices = list(range(83, 90))
        bad_indices_str = "83, 84, ..., 88, 89"
        with self.assertRaisesRegex(
            ValueError, rf"Selected index\(es\): \[{bad_indices_str}\] out of range"
        ):
            filter._indices_manager_to_perstructure(self.managers, bad_indices)

    def test_indices_perspecies(self):
        """Test the per-species, global-to-perstructure index transformation"""
        example_idces_perspecies = {1: [0, 1, 20], 8: [1, 2, 5], 6: [], 7: [0, 4]}
        # The order of the species should not matter;
        # the output is always sorted by species
        sps = [7, 6, 8, 1]
        perspecies_idces = filter._indices_perspecies_manager_to_perstructure(
            self.managers, example_idces_perspecies, sps
        )
        self.assertEqual(perspecies_idces, [[8, 9, 4], [], [9, 7, 0, 8], [], [0]])
        # Try it with an iterator (instead of a list) of species
        perspecies_idces = filter._indices_perspecies_manager_to_perstructure(
            self.managers, example_idces_perspecies, example_idces_perspecies.keys()
        )
        self.assertEqual(perspecies_idces, [[8, 9, 4], [], [9, 7, 0, 8], [], [0]])
        # And a set
        perspecies_idces = filter._indices_perspecies_manager_to_perstructure(
            self.managers, example_idces_perspecies, set(sps)
        )
        self.assertEqual(perspecies_idces, [[8, 9, 4], [], [9, 7, 0, 8], [], [0]])

    def test_indices_perspecies_out_of_range(self):
        """Test the per-species index transformer with out-of-range indices"""
        bad_idces_perspecies = {1: [], 8: [6], 7: [], 6: []}
        sps = [1, 6, 7, 8]
        with self.assertRaisesRegex(
            ValueError, r"Selected index\(es\): \[6\] for species 8 out of range"
        ):
            filter._indices_perspecies_manager_to_perstructure(
                self.managers, bad_idces_perspecies, sps
            )
        nonexistent_idces_perspecies = {1: [], 8: [], 7: [], 6: [], 12: [0]}
        sps = [1, 6, 7, 8, 12]
        with self.assertRaisesRegex(
            ValueError,
            r"Selected index\(es\): \[0\] for species 12 out of range "
            r"\(species does not appear to be present\)",
        ):
            filter._indices_perspecies_manager_to_perstructure(
                self.managers, nonexistent_idces_perspecies, sps
            )

    def test_indices_bad_species(self):
        """Test the global-to-perstructure selection with bad species lists"""
        # Missing species
        sps_bad = [1, 6, 8]
        sps_bad_re_str = r"\[1, 6, 8\]"
        with self.assertRaisesRegex(
            ValueError,
            f"^Atom of type 7 found but was not listed in sps: {sps_bad_re_str}$",
        ):
            filter._indices_perspecies_manager_to_perstructure(
                self.managers, {sp: [] for sp in sps_bad}, sps_bad
            )
        # Duplicated species
        sps_bad = [1, 6, 6, 7, 8]
        with self.assertRaisesRegex(
            ValueError,
            rf"^List of species contains duplicated entries: \[1, 6, 6, 7, 8\]",
        ):
            filter._indices_perspecies_manager_to_perstructure(
                self.managers, {sp: [0] for sp in sps_bad}, sps_bad
            )


class FilterTest:
    def abstractSetUp(self):
        example_frames = ase.io.read(
            "reference_data/inputs/small_molecules-20.json", ":"
        )
        global_species = set()
        for frame in example_frames:
            global_species.update([int(sp) for sp in frame.get_atomic_numbers()])
        repr = SOAP(
            max_radial=8,
            max_angular=6,
            interaction_cutoff=4.0,
            cutoff_smooth_width=1.0,
            gaussian_sigma_type="Constant",
            gaussian_sigma_constant=0.3,
            expansion_by_species_method="user defined",
            global_species=sorted(list(global_species)),
        )
        self.managers = repr.transform(example_frames)
        self.example_features = self.managers.get_features(repr)
        self.repr = repr
        self.example_frames = example_frames

    def test_sps_mode(self):
        """This test checks that `sample per species` mode returns the correct
        number of selected samples per each species
        """
        n_sparses = {1: 0, 6: 2, 7: 4, 8: 1}
        compressor = self._filter(self.repr, n_sparses, act_on="sample per species")

        x = compressor.select_and_filter(self.managers)

        # Checking that these parameters have not been set
        self.assertIsNone(compressor.selected_sample_ids)
        self.assertIsNone(compressor.selected_feature_ids_global)

        for sp in n_sparses:
            sp_idx = np.array(compressor.selected_sample_ids_by_sp[sp])
            self.assertEqual(len(sp_idx), n_sparses[sp])

        # Checking that the returned SparsePoints container has the right size
        n_sparse_total = sum(n_sparses[sp] for sp in n_sparses)
        self.assertEqual(x.size(), n_sparse_total)

    def test_sample_mode(self):
        """This test checks that `sample` mode returns the correct
        number of selected samples
        """
        n_sparses = self.example_features.shape[0] // 10
        compressor = self._filter(self.repr, n_sparses, act_on="sample")

        x = compressor.select_and_filter(self.managers)

        # Checking that these parameters have not been set
        self.assertIsNone(compressor.selected_feature_ids_global)
        self.assertIsNone(compressor.selected_sample_ids_by_sp)

        self.assertIsNotNone(compressor.selected_sample_ids)
        self.assertEqual(x.size(), n_sparses)

    def test_feature_mode(self):
        """This test checks that `feature` mode returns the correct
        number of selected feature
        """
        n_sparses = self.example_features.shape[1] // 100
        compressor = self._filter(self.repr, n_sparses, act_on="feature")
        compressor.select_and_filter(self.managers)

        self.assertIsNotNone(compressor.selected_feature_ids_global)
        idx = np.array(compressor.selected_feature_ids_global)
        self.assertEqual(len(idx), n_sparses)

        # Checking that these parameters have not been set
        self.assertIsNone(compressor.selected_sample_ids)
        self.assertIsNone(compressor.selected_sample_ids_by_sp)

    def test_missing_species(self):
        """This test checks that asking for a species which is not present
        raises an error.
        """
        n_sparses = {1: 0, 6: 2, 7: 4, 8: 1, 12: 3}
        compressor = self._filter(self.repr, n_sparses, act_on="sample per species")
        with self.assertRaisesRegex(
            ValueError,
            r"Found array with 0 sample\(s\) \(shape=\(0, 4480\)\) while a minimum of 1 is required.",
        ):
            compressor.select_and_filter(self.managers)

    def test_bad_mode(self):
        """This test checks that any mode other than `sample`,
        `sample per species`, and `feature` throws an error
        """
        with self.assertRaisesRegex(
            ValueError,
            '^"act_on" should be one of: "sample", "sample per species", or "feature"$',
        ):
            self._filter(self.repr, 1, act_on="bad mode")

    def test_new_n(self):
        """This test checks that `filter` can return an arbitrary number of
        selections strictly less than the number selected.
        """
        n_sparses = self.example_features.shape[1] // 100
        compressor = self._filter(self.repr, n_sparses, act_on="feature")
        compressor.select(self.managers)
        compressor.filter(self.managers, n_select=n_sparses // 10)

        with self.assertRaisesRegex(
            ValueError,
            rf"^It is only possible to filter {n_sparses} feature\(s\), you have requested {n_sparses * 10}$",
        ):
            compressor.filter(self.managers, n_select=n_sparses * 10)


class FPSTest(FilterTest, unittest.TestCase):
    def setUp(self):
        self._filter = FPSFilter
        self.abstractSetUp()


class CURTest(FilterTest, unittest.TestCase):
    def setUp(self):
        self._filter = CURFilter
        self.abstractSetUp()


if __name__ == "__main__":
    unittest.main(verbosity=2)
