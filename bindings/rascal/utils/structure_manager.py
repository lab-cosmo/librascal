import numpy as np

from ..lib import StructureManager, Adaptor


def ase2rascal(frame):
    """
    Convert ASE Atoms object to rascal's equivalent

    Parameters
    ----------
    frame : ase.Atoms
        Atomic structure

    Returns
    -------
    StructureManagerCenters
        base structure manager.
    """
    cell = frame.get_cell()
    positions = frame.get_positions()
    numbers = frame.get_atomic_numbers()
    pbc = frame.get_pbc().astype(int)

    managerC = StructureManager.Centers()
    managerC.update(np.array(positions.T, order='F'),
                    numbers.reshape(-1, 1),
                    np.array(cell.T, order='F'),
                    pbc.reshape(3, 1))
    return managerC


def get_strict_neighbourlist(frame, cutoff):
    """
    Generate a strict neighbourlist Structure Manager.

    Parameters
    ----------
    frame : ase.Atoms
        Atomic structure
    cutoff : float
        Cutoff radius for the strict neighbourlist

    Returns
    -------
    Strict_NeighbourList_Centers
        strict neighbourlist StructureManager
    """

    managerC = ase2rascal(frame)

    managerNL = Adaptor.NeighbourList_Centers(managerC, cutoff)
    managerNL.update()

    manager = Adaptor.Strict_NeighbourList_Centers(managerNL, cutoff)
    manager.update()

    return manager
