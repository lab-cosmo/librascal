from ..lib import sparsification

def fps(feature_matrix,Nselect,starting_index=None,flavour='simple'):
    """
    Farthest Point Sampling [1] routine using librascal.
    
    Parameters
    ----------
    feature_matrix : numpy.ndarray[float64[n, m], flags.c_contiguous]
        Feature matrix with n samples and m features.
    Nselect : int
        Number of sample to select
    starting_index : int, optional
        Index of the first sample to select (the default is None, which corresponds to starting_index == 0)
    flavour : str, optional
        Select which kind of FPS selection to perform:
        + 'simple' is the basic fps selection [1].
        + 'voronoi' is uses voronoi polyhedra to avoid some distance computations.
    
    Returns
    -------
    sparse_indices : numpy.ndarray[float64[Nselect]]
        Selected indices refering to the feature_matrix order.


    .. [1] Ceriotti, M., Tribello, G. A., & Parrinello, M. (2013). 
        Demonstrating the Transferability and the Descriptive Power of Sketch-Map. 
        Journal of Chemical Theory and Computation, 9(3), 1521–1532. 
        https://doi.org/10.1021/ct3010563
    """
    
    if starting_index is None:
        starting_index = 0

    if flavour == 'simple':
        sparse_indices,sparse_minmax_d2 = \
                sparsification.fps(feature_matrix,Nselect,starting_index)
        
    elif flavour == 'voronoi':
        sparse_indices, sparse_minmax_d2, voronoi_indices, voronoi_r2 = \
                sparsification.fps_voronoi(feature_matrix,Nselect,starting_index)
    
    else:
        raise Exception('Unknown flavour {}'.format(flavour))

    return sparse_indices,sparse_minmax_d2

