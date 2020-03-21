import numpy as np


def update_frame(frame, i_at, disp):
    """"Displace atom number i_at by displacement disp.

    Parameters
    ----------
    frame: ase.Atoms
    i_at: int < len(frame)
    disp: array of len 3

    Returns
    -------
    a displaced copy of frame

    """

    ff = frame.copy()
    ff.positions[i_at] += np.asarray(disp)
    # make sure atoms have not gone out of the unit cell
    ff.wrap(eps=1e-10)
    return ff

def compute_displaced_kernel(kernel, rep, frame, i_atom, disp, X_pseudo):
    """Compute the sparse kernel for the displaced version of frame (displacing
    atom i_atom by disp) with the set of sparse points X_pseudo

    Parameters
    ----------
    kernel: rascal.Kernel in sparse kernel mode
    rep: a rascal representation
    frame: ase.Atoms
    i_atom: int < len(frame)
    disp: array of len 3
    X_pseudo: rascal.PseudoPoints

    Returns
    -------

    """
    frame_d = update_frame(frame, i_atom, disp)
    managers = rep.transform([frame_d])
    KNM_d = kernel(managers, X_pseudo, grad=(False, False))
    return KNM_d

def compute_numerical_kernel_gradient(kernel, representation, frame, X_pseudo, eps=1e-5):
    """Compute the sparse kernel gradient with centered finite difference between the
    atomic structure frame and the sparse points X_pseudo.

    Parameters
    ----------
    kernel: rascal.Kernel in sparse kernel mode
    representation: a rascal representation
    frame: ase.Atoms
    X_pseudo: rascal.PseudoPoints
    eps=1e-5: finite difference displacement

    Returns
    -------
    np.array of dimension (len(frame)*3, X_pseudo.size())
    """
    disps = eps*np.array([[[1,0,0],[-1,0,0]],[[0,1,0],[0,-1,0]],[[0,0,1],[0,0,-1]]])
    KNM_der = np.zeros((len(frame)*3, X_pseudo.size()))

    for i_atom in range(len(frame)):
        for i_der, disp in enumerate(disps):
            KNM_p = compute_displaced_kernel(kernel, representation, frame, i_atom, disp[0], X_pseudo)
            KNM_m = compute_displaced_kernel(kernel, representation, frame, i_atom, disp[1], X_pseudo)
            KNM_der[i_atom*3+i_der] = ((KNM_p - KNM_m) / (2*eps)).sum(axis=0)
    return KNM_der

def compute_numerical_kernel_gradients(kernel, representation, frames, X_pseudo, eps=1e-5):
    """Compute the sparse kernel gradient with centered finite difference between the set of atomic structures frames and the sparse points X_pseudo.

    Parameters
    ----------
    kernel: rascal.Kernel in sparse kernel mode
    representation: a rascal representation
    frames: list of ase.Atoms
    X_pseudo: rascal.PseudoPoints
    eps=1e-5: finite difference displacement

    Returns
    -------
    np.array of dimension (len(frame)*3, X_pseudo.size())
    """
    KNM = []
    for frame in frames:
        KNM.append(compute_numerical_kernel_gradient(kernel, representation, frame, X_pseudo, eps))
    KNM = np.vstack(KNM)
    return KNM