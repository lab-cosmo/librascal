import inspect
from IPython.display import Markdown, display, clear_output
import ipywidgets as widgets
import numpy as np
import ase

def compute_representation(representation,frames):
    expansions = representation.transform(frames)
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
    model = KRR(zeta, weights, representation, features)
    return model,kernel


def show_readme_for_button(b):
    clear_output()
    display(Markdown(str(''.join(open('../README.rst')))))

def readme_button():
    button = widgets.Button(description="Shown README")
    output = widgets.Output()
    display(button, output)

    button.on_click(show_readme_for_button)

def make_buttons():
    import tutorial_utils as tutorial_utils

    funcs = {f[0]:f[1] for f in inspect.getmembers(tutorial_utils) if inspect.isfunction(f[1]) and 'button' not in f[0] and f[0] not in ['display','clear_output']}
    button = widgets.ToggleButtons(
            options = funcs.keys(),\
            description = "Show: ",\
            button_style='',\
            )

    def disp_func(a):
        clear_output()
        make_buttons()
        lines = inspect.getsource(funcs[button.value])
        display(Markdown('```python\n'+str(''.join(lines))+'```'))
    display(button)
    button.observe(disp_func, 'value')
