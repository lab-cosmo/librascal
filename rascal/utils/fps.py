from ..lib import sparsification
import numpy as np

def fps(feature_matrix, n_select, starting_index=None,
        method='simple', restart = None):
    """
    Farthest Point Sampling [1] routine using librascal.

    Parameters
    ----------
    feature_matrix : numpy.ndarray[float64[n, m], flags.c_contiguous]
        Feature matrix with n samples and m features.
    n_select : int
        Number of sample to select
    starting_index : int, optional
        Index of the first sample to select (the default is None,
                           which corresponds to starting_index == 0)
    restart : dictionary, optional (only valid for method="simple")
        the return value of a previous call to FPS. the selection
        will continue where it has been left off. the dictionary
        should contain at least "fps_indices"; if "fps_minmax_d2" and
        "fps_hausdorff_d2" are present, they will be used to restart,
        otherwise they'll be recomputed

    method : str, optional
        Select which kind of FPS selection to perform:
        + 'simple' is the basic fps selection [1].
        + 'voronoi' uses voronoi polyhedra to avoid some
                    distance computations.

    Returns
    -------
    A dictionary containing the following quantities
    "fps_indices" : numpy.ndarray[int[n_select]]
        Selected indices refering to the feature_matrix order.
    "fps_minmax_d2": numpy.ndarray[float64[n_select]]
        MIN-MAX distance^2 at each step in the selection
    "fps_hausforff_d2": numpy.ndarray[float64[n]]
        array of Hausdorff distances between the n points and the
        n_select FPS points



    .. [1] Ceriotti, M., Tribello, G. A., & Parrinello, M. (2013).
        Demonstrating the Transferability and the Descriptive Power of Sketch-Map.
        Journal of Chemical Theory and Computation, 9(3), 1521–1532.
        https://doi.org/10.1021/ct3010563
    """

    if starting_index is None:
        starting_index = 0


    return_dict = {}
    if method == 'simple':
        if restart is None:
            sparse_indices,sparse_minmax_d2,lmin_d2 = \
             sparsification.fps(feature_matrix,n_select,starting_index)
        else:
            res_tuple = ( restart["fps_indices"],
                ( restart["fps_minmax_d2"] if "fps_minmax_d2" in restart
                   else np.zeros(0, float) ),
                ( restart["fps_hausdorff_d2"] if "fps_hausdorff_d2" in restart
                   else np.zeros(0, float) ) )
            sparse_indices,sparse_minmax_d2,lmin_d2 = \
             sparsification.fps(feature_matrix,n_select,
              starting_index, res_tuple)

    elif method == 'voronoi':
        sparse_indices, sparse_minmax_d2, lmin_d2, \
            voronoi_indices, voronoi_r2 = \
            sparsification.fps_voronoi(feature_matrix,
                n_select,starting_index)
        return_dict["fps_voronoi_indices"] = voronoi_indices
        return_dict["fps_voronoi_r2"] = voronoi_r2

    else:
        raise Exception('Unknown FPS algorithm {}'.format(method))

    return_dict["fps_indices"] = sparse_indices
    return_dict["fps_minmax_d2"] = sparse_minmax_d2
    return_dict["fps_hausdorff_d2"] = lmin_d2

    return return_dict

