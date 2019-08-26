from rascal.utils import fps
import numpy as np

import unittest


def is_sorted(a):
    for i in range(a.size-1):
        if a[i+1] > a[i]:
            return False
    return True


class TestFPS(unittest.TestCase):

    def setUp(self):
        """
          Sets up a random feature matrix with 10000 points and
          dimensionality 100
        """
        np.random.seed(1234)
        self.x = np.random.normal(size=(10000, 100))

    def fps_return_match(self, r1, r2):
        """
          Check if return values match.
        """
        for key in ["fps_indices", "fps_minmax_d2", "fps_hausdorff_d2"]:
            self.assertTrue(np.allclose(r1[key], r2[key]))

    def test_simple_fps(self):
        """
          Just runs a small FPS run and see if the results make sense.
        """

        r_fps = fps(self.x, 100, 0, method="simple")
        self.assertTrue("fps_indices" in r_fps)
        self.assertTrue("fps_minmax_d2" in r_fps)
        self.assertTrue("fps_hausdorff_d2" in r_fps)
        self.assertTrue(r_fps["fps_hausdorff_d2"].max() <=
                        r_fps["fps_minmax_d2"][-2])
        self.assertTrue(is_sorted(r_fps["fps_minmax_d2"]))

    def test_voronoi_fps(self):
        """
          Just runs a small FPS run with voronoi method
          and see if the results make sense.
        """

        r_fps = fps(self.x, 100, 0, method="voronoi")
        self.assertTrue("fps_indices" in r_fps)
        self.assertTrue("fps_minmax_d2" in r_fps)
        self.assertTrue("fps_hausdorff_d2" in r_fps)
        self.assertTrue("fps_voronoi_indices" in r_fps)
        self.assertTrue("fps_voronoi_r2" in r_fps)
        self.assertTrue(r_fps["fps_hausdorff_d2"].max() <=
                        r_fps["fps_minmax_d2"][-2])
        self.assertTrue(is_sorted(r_fps["fps_minmax_d2"]))

    def test_fps_consistency(self):
        """
          Checks if different FPS methods give consistent results.
        """

        ref = fps(self.x, 100, 0, method="simple")
        r_fps = fps(self.x, 100, 0, method="voronoi")

        self.fps_return_match(r_fps, ref)

    def test_fps_restart(self):
        """
          Checks if FPS restart works as intended.
        """

        # This is the reference selection
        ref = fps(self.x, 100, 0, method="simple")

        # now we do this in two stages, by restarting
        # I first select 50
        r_fps = fps(self.x, 50, 0, method="simple")

        # I ask for 100 but provide the restart tuple.
        # Selection will continue from point 51
        r_fps = fps(self.x, 100, 0, method="simple", restart=r_fps)

        self.fps_return_match(r_fps, ref)

        # check if this works also when recomputing the distances
        r_fps = fps(self.x, 50, 0, method="simple")
        r_fps = fps(self.x, 100, 0, method="simple",
                    restart={"fps_indices": r_fps["fps_indices"]})
        self.fps_return_match(r_fps, ref)

        # You can also add new points (in this case we reuse
        # the same feature matrix) and continue the selection
        r_fps = fps(self.x, 50, 0, method="simple")
        xx = np.concatenate((self.x[r_fps["fps_indices"]], self.x))
        r_fps = fps(xx, 100, 0, method="simple",
                    restart={"fps_indices": np.asarray(range(50))})
        self. assertTrue(np.array_equal(r_fps["fps_indices"][:50],
                                        np.asarray(range(50))))
        self. assertTrue(np.array_equal(r_fps["fps_indices"][50:]-50,
                                        ref["fps_indices"][50:]))
