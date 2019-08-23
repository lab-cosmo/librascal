# Energy fitting

Approach used in Rascal can be used for calculating any of the rotationally invariant properties (such as formation energy), as well as rotationally covariant properties (such as dipole moment). This example is about fitting the formation energies of some small molecules. The dataset is located in the <a href="https://github.com/cosmo-epfl/librascal/tree/master/examples/data"><b>librascal/examples/data/</b></a> folder, the python notebook doing the calculation is 
<a href="https://github.com/cosmo-epfl/librascal/blob/master/examples/SOAP_example.ipynb"><b>librascal/examples/SOAP_example.ipynb </b></a>.

Let's take a look at the code and describe what it do step-by-step.

###Imports	
Just import of the necessary modules, including Rascal itself.
~~~~~~~~~~~~~{.py}
from matplotlib import pylab as plt
import os, sys
sys.path.insert(0,"../build/")
import time
import rascal
import json
import ase
from ase.io import read, writex
from ase.build import make_supercell
from ase.visualize import view
import numpy as np
from rascal.representations import SOAP
~~~~~~~~~~~~~

###Dataset
Here ASE is used to read the file with information about molecules, then the size of the feature matrix shown. The number of columns is the number of atoms in the set of molecules, the number of rows is depend on the max_radial and max_angular parameters. If the number of columns is big, sparsification needs to be done to reduce it.
~~~~~~~~~~~~~{.py}
frames = read('./data/small_molecules-1000.xyz',':')
hypers = dict(soap_type="PowerSpectrum",
              interaction_cutoff=3.5, 
              max_radial=6, 
              max_angular=6, 
              gaussian_sigma_constant=0.4,
              gaussian_sigma_type="Constant",
              cutoff_smooth_width=0.5,
              )
soap = SOAP(**hypers)
%time representation = soap.transform(frames)
X = representation.get_feature_matrix().T
X.shape
~~~~~~~~~~~~~

###Functions

Here we define the functions, which are needed for the further computations. Let's describe some of them. <br>
<b>compute_representation </b>- creates a feature manager class using information about structure<br>
<b>compute_kernel </b>- compute the matrix of global similarity using global cosine kernel with a power of zeta to itself<br>
<b>KRR.predict </b>- predicting the energy of the molecules in the set, using trained model<br>
<b>train_krr_model </b>- train the model basing on the given dataset and using global cosine kernel<br>

~~~~~~~~~~~~~{.py}
# Load the small molecules 
frames = read('./data/small_molecules-1000.xyz',':600')

def compute_representation(representation,frames):
    expansions = soap.transform(frames)
    return expansions

def compute_kernel(zeta, rep1, rep2=None):
    if rep2 is None:
        kernel = rep1.cosine_kernel_global(zeta)
    else:
        kernel = rep1.cosine_kernel_global(rep2,zeta)
    return kernel

def extract_energy(frames):
    prop = [[]]*len(frames)
    for ii,cc in enumerate(frames):
        prop[ii] = cc.info['dft_formation_energy_per_atom_in_eV']
    y = np.array(prop)
    return y

def split_dataset(frames, test_fraction, seed=10):
    N = len(frames)
    ids = np.arange(N)
    np.random.seed(seed)
    np.random.shuffle(ids)
    Ntrain = int(N*test_fraction)
    train = ids[:Ntrain]
    test = ids[Ntrain:]
    targets = extract_energy(frames)
    return [frames[ii] for ii in train],targets[train],[frames[ii] for ii in test],targets[test]

def get_mae(ypred,y):
    return np.mean(np.abs(ypred-y))
def get_rmse(ypred,y):
    return np.sqrt(np.mean((ypred-y)**2))
def get_sup(ypred,y):
    return np.amax(np.abs((ypred-y)))
def get_r2(y_pred,y_true):
    weight = 1
    sample_weight = None
    numerator = (weight * (y_true - y_pred) ** 2).sum(axis=0,dtype=np.float64)
    denominator = (weight * (y_true - np.average(
        y_true, axis=0, weights=sample_weight)) ** 2).sum(axis=0,dtype=np.float64)
    output_scores = 1 - (numerator / denominator)
    return np.mean(output_scores)

score_func = dict(
    MAE=get_mae,
    RMSE=get_rmse,
    SUP=get_sup,
    R2=get_r2,
)

def get_score(ypred,y):
    scores = {}
    for k,func in score_func.items():
        scores[k] = func(ypred,y)
    return scores

class KRR(object):
    def __init__(self,zeta,weights,representation,X):
        self.weights = weights
        self.representation = representation
        self.zeta = zeta
        self.X = X
        
    def predict(self,frames):
        features = compute_representation(self.representation,frames)
        kernel = compute_kernel(self.zeta , self.X, features)
        return np.dot(self.weights, kernel)
    
def train_krr_model(zeta,Lambda,representation,frames,y,jitter=1e-8):
    features = compute_representation(representation,frames)
    kernel = compute_kernel(zeta,features)
    # adjust the kernel so that it is properly scaled
    delta = np.std(y) / np.mean(kernel.diagonal())
    kernel[np.diag_indices_from(kernel)] += Lambda**2 / delta **2 + jitter
    # train the krr model
    weights = np.linalg.solve(kernel,y)
    model = KRR(zeta, weights,representation, features)
    return model,kernel
~~~~~~~~~~~~~

###Full spectrum

Here the full (with radial and angular parts) energies are computed. Let's describe the parameters of the soap descriptor, defined in the "hypers" dictionary.<br>
<b>interaction_cutoff</b> -  Maximum pairwise distance for atoms to be considered in expansion<br>
<b>max_radial</b> - number of radial basis functions<br>
<b>max_angular</b> - highest angular momentum number in the expansion<br>
<b>gaussian_sigma_constant</b> - specifies the atomic Gaussian widths, in the case where they're fixed.<br>
<b>gaussian_sigma_type</b> - how the Gaussian atom sigmas (smearing widths) are allowed to vary -- fixed ('Constant'), by species ('PerSpecies'), or by distance from the central atom ('Radial')<br>
<b>cutoff_smooth_width</b> - the distance over which the the interaction is smoothed to zero<br>

~~~~~~~~~~~~~{.py}
hypers = dict(soap_type="PowerSpectrum",
              interaction_cutoff=3.5, 
              max_radial=6, 
              max_angular=6, 
              gaussian_sigma_constant=0.4,
              gaussian_sigma_type="Constant",
              cutoff_smooth_width=0.5,
              )
soap = SOAP(**hypers)

frames_train, y_train, frames_test, y_test = split_dataset(frames,0.8)

zeta = 2
Lambda = 5e-3
krr,k = train_krr_model(zeta, Lambda, soap, frames_train, y_train)

y_pred = krr.predict(frames_test)
get_score(y_pred, y_test)

plt.scatter(y_pred, y_test, s=3)
plt.axis('scaled')
plt.xlabel('DFT energy / (eV/atom)')
plt.ylabel('Predicted energy / (eV/atom)')
~~~~~~~~~~~~~
The result of this block is:
![Prediction of the energy and the real energy of each molecule](../../src/R1s.png)
The result is quite good. One can try to change the train dataset to see how it affects the precision of the result. 

###Radial spectrum

Here we compute the energy, supposing the angular component to be zero.
~~~~~~~~~~~~~{.py}
hypers = dict(soap_type="RadialSpectrum",
              interaction_cutoff=3.5, 
              max_radial=6, 
              max_angular=0, 
              gaussian_sigma_constant=0.4,
              gaussian_sigma_type="Constant",
              cutoff_smooth_width=0.5,
              )
soap = SOAP(**hypers)

frames_train, y_train, frames_test, y_test = split_dataset(frames,0.8)

zeta = 2
Lambda = 5e-4
krr,k = train_krr_model(zeta, Lambda, soap, frames_train, y_train)

y_pred = krr.predict(frames_test)
get_score(y_pred, y_test)

plt.scatter(y_pred, y_test, s=3)
plt.axis('scaled')
plt.xlabel('DFT energy / (eV/atom)')
plt.ylabel('Predicted energy / (eV/atom)')
~~~~~~~~~~~~~
Comparison of full and radial spectrum:
![Comparison of full and radial spectrum](../../src/Comps.png)

It can be seen that the two spectres are quite similar, but the radial spectrum is much more simple to compute (as feature matrix is much smaller and the set of spherical harmonics doesn't have to be computed). It is quite an inteseting fact, but, unfortunately, this feature is probably not generalizable and should be just the feature of this particular dataset.

###Map of the dataset

Here we use sklearn to do <a href="https://en.wikipedia.org/wiki/Kernel_principal_component_analysis"><b>kernel principal component analysis</b></a>. 
~~~~~~~~~~~~~{.py}
def compute_representation(representation,frames):
    expansions = soap.transform(frames)
    return expansions

def compute_kernel(zeta, rep1, rep2=None):
    if rep2 is None:
        kernel = rep1.cosine_kernel_global(zeta)
    else:
        kernel = rep1.cosine_kernel_global(rep2,zeta)
    return kernel

def link_ngl_wdgt_to_ax_pos(ax, pos, ngl_widget):
    from matplotlib.widgets import AxesWidget
    from scipy.spatial import cKDTree
    r"""
    Initial idea for this function comes from @arose, the rest is @gph82 and @clonker
    """
    
    kdtree = cKDTree(pos)
    #assert ngl_widget.trajectory_0.n_frames == pos.shape[0]
    x, y = pos.T
    
    lineh = ax.axhline(ax.get_ybound()[0], c="black", ls='--')
    linev = ax.axvline(ax.get_xbound()[0], c="black", ls='--')
    dot, = ax.plot(pos[0,0],pos[0,1], 'o', c='red', ms=7)

    ngl_widget.isClick = False
    
    def onclick(event):
        linev.set_xdata((event.xdata, event.xdata))
        lineh.set_ydata((event.ydata, event.ydata))
        data = [event.xdata, event.ydata]
        _, index = kdtree.query(x=data, k=1)
        dot.set_xdata((x[index]))
        dot.set_ydata((y[index]))
        ngl_widget.isClick = True
        ngl_widget.frame = index
    
    def my_observer(change):
        r"""Here comes the code that you want to execute
        """
        ngl_widget.isClick = False
        _idx = change["new"]
        try:
            dot.set_xdata((x[_idx]))
            dot.set_ydata((y[_idx]))            
        except IndexError as e:
            dot.set_xdata((x[0]))
            dot.set_ydata((y[0]))
            print("caught index error with index %s (new=%s, old=%s)" % (_idx, change["new"], change["old"]))
    
    # Connect axes to widget
    axes_widget = AxesWidget(ax)
    axes_widget.connect_event('button_release_event', onclick)
    
    # Connect widget to axes
    ngl_widget.observe(my_observer, "frame", "change")

# Load the small molecules 
frames = read('./data/small_molecules-1000.xyz',':600')
hypers = dict(soap_type="PowerSpectrum",
              interaction_cutoff=3.5, 
              max_radial=6, 
              max_angular=6, 
              gaussian_sigma_constant=0.4,
              gaussian_sigma_type="Constant",
              cutoff_smooth_width=0.5,
              )
soap = SOAP(**hypers)
zeta = 2
features = compute_representation(soap, frames)
kernel = compute_kernel(zeta,features)
from sklearn.decomposition import KernelPCA
kpca = KernelPCA(n_components=2,kernel='precomputed')
kpca.fit(kernel)
X = kpca.transform(kernel)
plt.scatter(X[:,0],X[:,1],s=3)
~~~~~~~~~~~~~

The result of this block is:
![KernalPCA](../../src/PCAs.png)
It shows how the structures is located in the abstract 2D map, where similar structures are located near to each other, and the very different ones far from each other. 
