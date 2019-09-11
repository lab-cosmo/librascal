import inspect
import os
from IPython.display import Markdown, display, clear_output
import ipywidgets as widgets
import numpy as np
import ase
from ase.io import read
try:
    from rascal.representations import SphericalInvariants as SOAP
except:
    from rascal.representations import SOAP
from matplotlib import pyplot as plt
import time
from time import strftime
from tqdm import tqdm_notebook as tqdm

style = {'description_width': 'initial'}
hyper_vals = {"soap_type": ["PowerSpectrum", "RadialSpectrum"],
            "interaction_cutoff":[2.0,5.0],
            "max_radial": [0, 10],
            "max_angular": [0, 10],
            "gaussian_sigma_constant": [0.0,1.0],
            "gaussian_sigma_type":["Constant"],
            "cutoff_smooth_width":[0.0,1.0],
            }
hyper_dict = {"Power Spectrum":dict(soap_type="PowerSpectrum",
                              interaction_cutoff=3.5,
                              max_radial=2,
                              max_angular=1,
                              gaussian_sigma_constant=0.5,
                              gaussian_sigma_type="Constant",
                              cutoff_smooth_width=0.,
                              ),\
          "Full Power Spectrum":dict(soap_type="PowerSpectrum",
                              interaction_cutoff=3.5,
                              max_radial=6,
                              max_angular=6,
                              gaussian_sigma_constant=0.4,
                              gaussian_sigma_type="Constant",
                              cutoff_smooth_width=0.5,
                              ),
          "Radial Spectrum": dict(soap_type="RadialSpectrum",
                                  interaction_cutoff=3.5,
                                  max_radial=6,
                                  max_angular=0,
                                  gaussian_sigma_constant=0.4,
                                  gaussian_sigma_type="Constant",
                                  cutoff_smooth_width=0.5,
                                  )
          }
known_properties = dict(ENERGY="atomic", dft_formation_energy_per_atom_in_eV="global")
def split_dataset(frames, test_fraction, seed=10):
    N = len(frames)
    ids = np.arange(N)
    np.random.seed(seed)
    np.random.shuffle(ids)
    Ntrain = int(N*test_fraction)
    train = ids[:Ntrain]
    test = ids[Ntrain:]
    targets = extract_energy(frames)
    return np.array(train), np.array(test)#[frames[ii] for ii in train],targets[train],[frames[ii] for ii in test],targets[test]
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
def get_score(ypred,y):
    score_func = dict(
        MAE=get_mae,
        RMSE=get_rmse,
        SUP=get_sup,
        R2=get_r2,
    )
    scores = {}
    for k,func in score_func.items():
        scores[k] = func(ypred,y)
    return scores
def compute_representation(representation,frames):
    expansions = representation.transform(frames)
    return expansions
def compute_kernel(zeta, rep1, rep2=None, kernel_type = 'global'):
    if(kernel_type=='atomic'):
        kernel_function = cosine_kernel_atomic
    else:
        kernel_function = cosine_kernal_global
    if rep2 is not None:
        return rep1.kernel_function(rep2, zeta)
    else:
        return rep1.kernel_function(zeta)
class KRR(object):
    def __init__(self,zeta,weights,representation,X, kernel_type):
        self.weights = weights
        self.representation = representation
        self.zeta = zeta
        self.X = X
        self.kernel_type=kernel_type

    def predict(self,frames):
        features = compute_representation(self.representation,frames)
        kernel = compute_kernel(self.zeta, self.X, features, self.kernel_type)
        return np.dot(self.weights, kernel)
def extract_property(frames, property='energy'):
    try:
        return np.array([cc.info[property] for cc in frames])
    except:
        print(frames[0].info)
        raise KeyError("{} is not a property in the given frames".format(property))
def readme_button():
    def show_readme_for_button(b):
        clear_output()
        display(Markdown(str(''.join(open('../README.rst')))))
    button = widgets.Button(description="Show README", style=style)
    output = widgets.Output()
    display(button, output)

    button.on_click(show_readme_for_button)
def make_util_buttons():
    funcs = {'compute_kernel': compute_kernel,
             'compute_representation': compute_representation,
             'extract_energy': extract_energy,
             'get_mae': get_mae,
             'get_r2': get_r2,
             'get_rmse': get_rmse,
             'get_score': get_score,
             'get_sup': get_sup,
             'split_dataset': split_dataset,
             'train_krr_model': train_krr_model
            }

    def disp_func(a):
        clear_output()
        make_util_buttons()
        lines = inspect.getsource(funcs[button.value])
        display(Markdown('```python\n'+str(''.join(lines))+'```'))

    button = _button_template_(funcs.keys(), "Show: ", disp_func)
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
def _button_template_(options, description, disp_func, default=None):
    button = widgets.ToggleButtons(
        options = options,\
        description = description,
        button_style='',\
        style=style,\
        value = options[0] if default==None else default
        )
    display(button)
    button.observe(disp_func, 'value')

    return button

# def setup_SOAP():
class SOAP_tutorial(object):
    def __init__(self, input_file='small_molecules-1000.xyz',perc_test=0.8, interactive = False, verbose=True, hyperparameters = dict(**hyper_dict['Power Spectrum']), number_of_frames=None, property=None):
        self.perc_test = perc_test
        self.zeta = 2
        self.Lambda = 5e-3

        self.hyperparameters = hyperparameters
        self.verbosity_wrap = lambda s: (print('\r') if not verbose else print(s))

        file_options = [f for f in os.listdir('./data') if f.endswith('xyz')]
        if(input_file!=None):
            # self.verbosity_wrap("Our input file is {}.\n".format(input_file))
            file_options.insert(0, input_file)
            if(input_file in file_options[1:]):
                file_options.pop(file_options[1:].index(input_file)+1)
            # self.verbosity_wrap("Our hyperparameters are \n\t{}".format('\n\t'.join(['{:<30} {:<10}'.format(k, self.hyperparameters[k]) for k in self.hyperparameters])))
        self.input_file = file_options[0] if input_file==None else input_file
        self.frames=np.array(read('./data/'+self.input_file, ":"))
        self.train_idx, self.test_idx = None, None

        self.sliders = {val:
                        widgets.FloatSlider(value = self.hyperparameters.get(val, hyper_vals[val][0]),
                                    min = hyper_vals[val][0],
                                    max = hyper_vals[val][1],
                                    description = val.replace("_"," ").title(),
                                    continuous_update = True,
                                    step = (hyper_vals[val][1]-hyper_vals[val][0])/20.,\
                                    style=style,\
                                    )
                        if isinstance(hyper_vals[val][0],float) else
                        widgets.IntSlider(value = self.hyperparameters.get(val, hyper_vals[val][0]),
                                    min = hyper_vals[val][0],
                                    max = hyper_vals[val][1],
                                    description = val.replace("_"," ").title(),
                                    continuous_update = True,
                                    step = 1,
                                    style=style,\
                                    )
                        if isinstance(hyper_vals[val][0],int) and hyper_vals[val][0]!=True else
                        widgets.Dropdown(options=hyper_vals[val],
                                         style=style,\
                                         value = self.hyperparameters.get(val, hyper_vals[val][0]),
                                         description=val.replace("_"," ").title())
                        for val in hyper_vals}

        self.input_button = _button_template_(file_options, "Input File: ", self.get_input)

        self.preset_button = _button_template_(list(hyper_dict.keys()), "SOAP Presets: ", self.disp_func)

        self.properties = {prop: extract_property(self.frames, prop) for prop in self.frames[0].info if prop in known_properties}
        self.sliders['number_of_frames'] = widgets.IntSlider(value=len(self.frames) if number_of_frames==None else number_of_frames, min = 0, max = len(self.frames), \
                                                    description="Number of Frames", step=1, style=style)
        self.sliders['property'] = widgets.Dropdown(value=list(self.properties.keys())[0] if property==None else property, options=list(self.properties.keys()), \
                                                    description="Property to ML", style=style)
        if(interactive):
            for s in self.sliders:
                display(self.sliders[s])
                self.sliders[s].observe(lambda change: self.change_func(change), names='value')
        self.krr = {prop:None for prop in self.properties}
        self.trained={prop:False for prop in self.properties}
    def reset_ML(self, inp_change=False):
        if(inp_change):
            self.properties = {prop: extract_property(self.frames, prop) for prop in self.frames[0].info if  prop in known_properties}
            self.sliders['number_of_frames'].max = len(self.frames)
            self.sliders['number_of_frames'].value = len(self.frames)
            self.sliders['property'].value = list(self.properties.keys())[0]
            self.sliders['property'].options=list(self.properties.keys())
        self.krr = {prop:None for prop in self.properties}
        self.trained={prop:False for prop in self.properties}
    def change_func(self,change):
        slider = self.sliders[change['owner'].description.replace(' ','_').lower()]
        if(slider.description.replace(' ','_').lower() in self.hyperparameters):
            self.hyperparameters[slider.description.replace(' ','_').lower()] = change['new']
        slider.value = change['new']
        self.reset_ML()
    def disp_func(self,a):
        self.hyperparameters = hyper_dict[self.preset_button.value]
        for s in self.sliders:
            if(s in self.hyperparameters):
                self.sliders[s].value = self.hyperparameters[s]
        self.reset_ML()
    def get_input(self,a):
        self.input_file = self.input_button.value
        self.frames = np.array(read('./data/'+ self.input_file, ":"))
        self.reset_ML(True)
    def compute_SOAP(self, frames = None):
        if(frames==None):
            self.train_idx, self.test_idx = split_dataset(self.frames[:self.sliders['number_of_frames'].value],self.perc_test)
            frames=self.frames[self.train_idx]

        self.verbosity_wrap("\nWe'll now construct the SOAP representation of the first frame by inputting our hyperparameters into rascal.representations.SphericalInvariants and use this to transform our frames.")
        soap = SOAP(**self.hyperparameters)

        compute_representation(soap, frames)
    def train_krr_model(self, jitter=1e-8):
        representation = SOAP(**self.hyperparameters)

        self.train_idx, self.test_idx = split_dataset(self.frames[:self.sliders['number_of_frames'].value],self.perc_test)
        self.verbosity_wrap("First, I am going to separate my dataset:\n\tTraining Set: {} pts ({}%)\n\tTesting Set: {} pts ({}%)".format(
                                                                  len(self.train_idx), \
                                                                  round(100*(self.perc_test)), \
                                                                  len(self.test_idx), \
                                                                  round(100*(1-self.perc_test))))
        t = time.time()

        self.verbosity_wrap("\nNow we will compute the SOAP representation of our training frames.")

        features = compute_representation(representation, self.frames[self.train_idx])

        self.verbosity_wrap('In this run, computing the SOAP vectors took {} seconds/frame'.format(round((time.time()-t)/len(self.train_idx),8)))
        self.verbosity_wrap("\nNext we find the kernel for our training model.\n(This step may take a few minutes for larger training sets.)")

        time.sleep(.5)
        kernel = compute_kernel(self.zeta, features, kernel_type=known_properties[self.sliders['property'].value])
        delta = np.std(self.properties[self.sliders['property'].value][self.train_idx]) / np.mean(kernel.diagonal())
        kernel[np.diag_indices_from(kernel)] += self.Lambda**2 / delta **2 + jitter

        self.verbosity_wrap("\nWe will adjust the diagonals of our kernel by {} so that it is properly scaled.".format(self.Lambda**2 / delta **2 + jitter))
        self.verbosity_wrap("\nNow we can take this kernel to compute the weights of our KRR.")

        weights = np.linalg.solve(kernel,self.properties[self.sliders['property'].value][self.train_idx])
        model = KRR(self.zeta, weights, representation, features, kernel_type=known_properties[self.sliders['property'].value])
        self.krr[self.sliders['property'].value], k = model, kernel
        self.trained[self.sliders['property'].value] = True
        return self.krr[self.sliders['property'].value], k

    def plot_prediction(self):
        if(self.trained[self.sliders['property'].value]==False):
            self.train_model(self.sliders['property'].value)
        y_pred = self.krr[self.sliders['property'].value].predict(self.frames[self.test_idx])
        print(get_score(y_pred, self.properties[self.sliders['property'].value][self.test_idx]))
        plt.scatter(y_pred, self.properties[self.sliders['property'].value][self.test_idx], s=3)
        plt.axis('scaled')
        plt.xlabel('DFT energy / (eV/atom)')
        plt.ylabel('Predicted energy / (eV/atom)')
        plt.gca().set_aspect('equal')
        plt.show()


if __name__=="__main__":
    make_util_buttons()