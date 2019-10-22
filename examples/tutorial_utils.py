import os
import inspect
from IPython.display import Markdown, display, clear_output
import ipywidgets as widgets
import numpy as np
import ase
from ase.io import read
try:
    from rascal.representations import SphericalInvariants as SOAP
except:
    from rascal.representations import SOAP
from rascal.models import Kernel
from matplotlib import pyplot as plt
import time
from time import strftime
from tqdm import tqdm_notebook as tqdm
from threading import Thread, Timer
import urllib.request
from queue import Queue
import multiprocessing

style = {'description_width': 'initial'}
hyper_vals = {
# "soap_type": dict(options = ["PowerSpectrum", "RadialSpectrum"], name = "Body Order"),
            "interaction_cutoff":dict(options = [2.0,5.0],name = r"$r_{cut}$"),
            "max_radial":dict(options =  [0, 10],name = r"$n_{max}$"),
            "max_angular": dict(options = [0, 10],name = r"$l_{max}$"),
            "gaussian_sigma_constant": dict(options = [0.01,1.0],name = r"$\sigma$"),
            "gaussian_sigma_type":dict(fixed="Constant"),
            "cutoff_smooth_width":dict(fixed=0.5),
            }
hyper_dict = {"Minimal Power Spectrum":dict(
                              # soap_type="PowerSpectrum",
                              interaction_cutoff=3.5,
                              max_radial=2,
                              max_angular=1,
                              gaussian_sigma_constant=0.5,
                              gaussian_sigma_type="Constant",
                              cutoff_smooth_width=0.,
                              ),\
          "Power Spectrum":dict(
                              # soap_type="PowerSpectrum",
                              interaction_cutoff=3.5,
                              max_radial=6,
                              max_angular=6,
                              gaussian_sigma_constant=0.4,
                              gaussian_sigma_type="Constant",
                              cutoff_smooth_width=0.5,
                              ),
          "Radial Spectrum": dict(
                                  # soap_type="RadialSpectrum",
                                  interaction_cutoff=3.5,
                                  max_radial=6,
                                  max_angular=0,
                                  gaussian_sigma_constant=0.4,
                                  gaussian_sigma_type="Constant",
                                  cutoff_smooth_width=0.5,
                                  )
          }
known_properties = dict(CS="atomic", dft_formation_energy_per_atom_in_eV="global")
ignore= ['Natoms', 'numbers','Name','positions','cutoff','nneightol',"NAME"]

def markdown_table_from_dict(d, headers):
    return '<table>  <thead><tr>{}</tr></thead><tbody>'.format(''.join(['<th>{}</th>'.format(h) for h in headers]))+''.join(['<tr><td>{}</td>{}</tr>'.format(key, ''.join(['<td>{}</td>'.format(v) for v in d[key]])) for key in d])+'</tbody></table>'
def split_dataset(frames, train_fraction, seed=10):
    N = len(frames)
    ids = np.arange(N)
    np.random.seed(seed)
    np.random.shuffle(ids)
    Ntrain = int(N*train_fraction)
    train = ids[:Ntrain]
    test = ids[Ntrain:]
    return np.array(train), np.array(test)#[frames[ii] for ii in train],targets[train],[frames[ii] for ii in test],targets[test]
def get_r2(y_pred,y_true):
    numerator = ((y_true - y_pred) ** 2).sum(axis=0,dtype=np.float64)
    denominator = ((y_true - np.average(y_true, axis=0)) ** 2).sum(axis=0,dtype=np.float64)
    output_scores = 1 - numerator / denominator
    return np.mean(output_scores)
def get_score(ypred,y):
    return {
                "SUP": [np.amax(np.abs((ypred-y)))],
                "MAE": [np.mean(np.abs(ypred-y))],
                "RMSD":[np.sqrt(np.mean((ypred-y)**2))],
                r"$R^2$" : [get_r2(ypred, y)]
                }
def compute_kernel(representation, zeta, feature1, feature2=None, kernel_type = 'global', time_estimate=0):
    if(kernel_type=='atomic'):
        kernel = Kernel(representation, name='Cosine',target_type='Atom',zeta=zeta)
    else:
        kernel = Kernel(representation, name='Cosine',target_type='Structure',zeta=zeta)
    # if(time_estimate!=0):
    #     manager = multiprocessing.Manager()
    #     return_dict = manager.dict()
    #     # timer = Thread(target=tqdm_timer, args=(time_estimate,))
    #     timer = multiprocessing.Process(target=tqdm_timer, args=(time_estimate,))
    #     timer.start()
    #     # kthread = Thread(target=kernel_function, args=(rep2, zeta,) if rep2!=None else (zeta,))
    #     kthread = multiprocessing.Process(target=kernel_function, args=(rep2, zeta,) if rep2!=None else (zeta,))
    #     kthread.start()
    #     kthread.join()
    #     print(return_dict.values())
    #     return return_dict.values()
    # else:
    if feature2 is not None:
        return kernel(feature2)
    else:
        return kernel(feature1)

class KRR(object):
    def __init__(self,zeta,weights,representation,X, kernel_type):
        self.weights = weights
        self.representation = representation
        self.zeta = zeta
        self.X = X
        self.kernel_type=kernel_type

    def predict(self,frames):
        features = self.representation.transform(frames)
        kernel = compute_kernel(self.representation, self.zeta, self.X, features, self.kernel_type)
        return np.dot(self.weights, kernel)
def extract_property(frames, property='energy'):
    if(property in frames[0].info):
        return np.array([cc.info[property] for cc in frames])
    elif(property in frames[0].arrays):
        return np.array([cc.arrays[property] for cc in frames])
    else:
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
def tqdm_timer(seconds):
    for i in tqdm(range(seconds), total=seconds):
        time.sleep(1)
class SOAP_tutorial(object):

    csd500_filename = "CSD-500.xyz"
    csd500_download_link = (
            "https://static-content.springer.com/esm/" +
            "art%3A10.1038%2Fs41467-018-06972-x/MediaObjects/" +
            "41467_2018_6972_MOESM4_ESM.txt")

    def __init__(self,
                 input_file='./reference_data/small_molecules-1000.xyz',
                 auto_download_csd500=True,
                 training_percentage=0.8,
                 interactive = False, verbose=True, hyperparameters=dict(
                     **hyper_dict['Power Spectrum']),
                 number_of_frames=None, property=None):
        self.zeta = 2
        self.Lambda = 5e-3

        self.hyperparameters = hyperparameters
        self.verbose = verbose
        self.verbosity_wrap = lambda s: (None if not verbose else display(Markdown(s)))

        if auto_download_csd500 and not os.path.isfile(os.path.join(
                "data", self.csd500_filename)):
            try:
                urllib.request.urlretrieve(
                        self.csd500_download_link,
                        os.path.join("data", self.csd500_filename))
            except urllib.error.URLError as e:
                sys.stderr.write("Couldn't download the CSD-500 dataset.\n")
                if hasattr(e, "reason"):
                    sys.stderr.write("Failed to reach the server. " +
                                     "Try checking your network connection?\n")
                    sys.stderr.write("Reason: " + e.reason + "\n")
                elif hasattr(e, "code"):
                    sys.stderr.write("HTTP error: {:d}\n".format(e.code))
        file_options = ['./reference_data/{}'.format(f) for f in os.listdir('./data') if f.endswith('xyz')]
        if(input_file!=None):
            file_options.insert(0, input_file)
            if(input_file in file_options[1:]):
                file_options.pop(file_options[1:].index(input_file)+1)
        self.input_file = file_options[0] if input_file==None else input_file
        self.frames=np.array(read(self.input_file, ":"))
        self.train_idx, self.test_idx = None, None

        self.sliders = {val:
                        widgets.FloatSlider(value = self.hyperparameters.get(val, hyper_vals[val]['options'][0]),
                                    min = hyper_vals[val]['options'][0],
                                    max = hyper_vals[val]['options'][1],
                                    description = hyper_vals[val]['name'],
                                    continuous_update = True,
                                    step = (hyper_vals[val]['options'][1]-hyper_vals[val]['options'][0])/20.,\
                                    style=style,\
                                    )
                        if isinstance(hyper_vals[val]['options'][0],float) else
                        widgets.IntSlider(value = self.hyperparameters.get(val, hyper_vals[val]['options'][0]),
                                    min = hyper_vals[val]['options'][0],
                                    max = hyper_vals[val]['options'][1],
                                    description = hyper_vals[val]['name'],
                                    continuous_update = True,
                                    step = 1,
                                    style=style,\
                                    )
                        if isinstance(hyper_vals[val]['options'][0],int) and hyper_vals[val]['options'][0]!=True else
                        widgets.Dropdown(options=hyper_vals[val]['options'],
                                         style=style,\
                                         value = self.hyperparameters.get(val, hyper_vals[val]['options'][0]),
                                         description=hyper_vals[val]['name'],
                                         )
                        for val in hyper_vals if 'fixed' not in hyper_vals[val]}


        self.properties = {prop: extract_property(self.frames, prop) for prop in sorted([*list(self.frames[0].info.keys()), *list(self.frames[0].arrays.keys())], key=lambda p: -int(p in known_properties)) if prop not in ignore}
        self.sliders['number_of_frames'] = widgets.IntSlider(value=int(len(self.frames)*0.2) if number_of_frames==None else number_of_frames, min = 1, max = len(self.frames), \
                                                    description="Number of Frames", step=1, style=style)
        self.sliders['property_to_ml'] = widgets.Dropdown(value=list(self.properties.keys())[0] if property==None else property, options=list(self.properties.keys()), \
                                                    description="Property to ML", style=style)
        self.sliders['kernel_type'] = widgets.Dropdown(value='atomic' if self.sliders['property_to_ml'].value not in known_properties else known_properties[self.sliders['property_to_ml'].value],
                                                    options=['atomic', 'global'] if self.sliders['property_to_ml'].value not in known_properties else [known_properties[self.sliders['property_to_ml'].value]], \
                                                    description="Kernel Type", style=style)
        self.sliders['training_percentage'] = widgets.FloatSlider(value = training_percentage,
                                                                    min = 0,
                                                                    max = 1,
                                                                    description = "Training Percentage",
                                                                    continuous_update = True,
                                                                    step = 0.05,
                                                                    style=style,\
                                                                    )
        if(interactive):
            self.input_button = widgets.Dropdown(options=[*file_options],#, "Other"],
                                         style=style,\
                                         value = self.input_file,
                                         description="Input File: ",
                                         )
            self.input_button.observe(self.get_input, names='value')
            display(self.input_button)
            self.preset_button = _button_template_(list(hyper_dict.keys()), "SOAP Presets: ", self.disp_func)

            slider_order = ['property_to_ml', 'kernel_type', 'number_of_frames', 'training_percentage',*list(hyper_vals.keys())]
            for s in slider_order:
                if(s in self.sliders):
                    display(self.sliders[s])
                    self.sliders[s].observe(lambda change: self.change_func(change), names='value')
        self.krr = {prop:None for prop in self.properties}
        self.trained={prop:False for prop in self.properties}
        self.est_frames, self.est_times = [], []
        self.pred_frames, self.pred_times = [], []
    def reset_ML(self, inp_change=False):
        if(inp_change):
            self.properties = {prop: extract_property(self.frames, prop) for prop in sorted([*list(self.frames[0].info.keys()), *list(self.frames[0].arrays.keys())], key=lambda p: -int(p in known_properties)) if prop not in ignore}
            self.sliders['number_of_frames'].max = len(self.frames)
            self.sliders['number_of_frames'].value = int(len(self.frames)*0.2)
            self.sliders['property_to_ml'].options=list(self.properties.keys())
            self.sliders['property_to_ml'].value = list(self.properties.keys())[0]

            self.sliders['kernel_type'].options  =['atomic', 'global'] if self.sliders['property_to_ml'].value not in known_properties else [known_properties[self.sliders['property_to_ml'].value]]
            self.sliders['kernel_type'].value = 'atomic' if self.sliders['property_to_ml'].value not in known_properties else known_properties[self.sliders['property_to_ml'].value]

        self.krr = {prop:None for prop in self.properties}
        self.trained={prop:False for prop in self.properties}
    def change_func(self,change):
        change['owner'].value = change['new']
        for s in self.hyperparameters:
            if(s in self.sliders):
                if(self.hyperparameters[s]!=self.sliders[s].value):
                    self.est_frames, self.est_times = [], []
                    self.hyperparameters[s] = self.sliders[s].value
        if(change['owner']==self.sliders['property_to_ml']):
            self.sliders['kernel_type'].options  =['atomic', 'global'] if self.sliders['property_to_ml'].value not in known_properties else [known_properties[self.sliders['property_to_ml'].value]]
            self.sliders['kernel_type'].value = 'atomic' if self.sliders['property_to_ml'].value not in known_properties else known_properties[self.sliders['property_to_ml'].value]

        self.reset_ML()
    def disp_func(self,a):
        for s in self.sliders:
            if(s in self.hyperparameters):
                self.hyperparameters[s] = hyper_dict[self.preset_button.value][s]
                self.sliders[s].value = self.hyperparameters[s]
        self.reset_ML()
        # timer=Thread(target=tqdm_timer, args=(e
    def get_input(self,a):
        inp = self.input_button.value# if self.input_button.value!='Other' else input("Which file do you want to load?")
        self.input_file = inp
        self.frames = np.array(read(self.input_file, ":"))
        self.reset_ML(True)
    def train_krr_model(self, jitter=1e-8):

        self.verbosity_wrap(markdown_table_from_dict({"Training Set": [int(self.sliders['number_of_frames'].value*(self.sliders['training_percentage'].value)),round(100*(self.sliders['training_percentage'].value))],
                                                      "Testing Set":[int(self.sliders['number_of_frames'].value*(1-self.sliders['training_percentage'].value)),round(100*(1-self.sliders['training_percentage'].value))]}, headers=["Partition","Number of Frames", "Percentage"]))

        self.train_idx, self.test_idx = split_dataset(self.frames[:self.sliders['number_of_frames'].value],self.sliders['training_percentage'].value)
        self.train_krr_model_func(self.train_idx, jitter=jitter)
    def train_krr_model_func(self, frame_idx, jitter=1e-8, pretend=False):
        if(pretend==False):
            self.output_params()
            verbosity_wrap = lambda s: self.verbosity_wrap(s)
        else:
            verbosity_wrap = lambda s: None

        verbosity_wrap("<br/>We will now train a model on {}.".format(self.sliders['property_to_ml'].value))

        representation = SOAP(**self.hyperparameters)

        training_properties = np.concatenate(self.properties[self.sliders['property_to_ml'].value][frame_idx])[:,0] if known_properties[self.sliders['property_to_ml'].value]=='atomic' else self.properties[self.sliders['property_to_ml'].value][frame_idx]

        verbosity_wrap("First, I am going to separate my dataset:")
        verbosity_wrap("<br/>Now we will compute the SOAP representation of our training frames.")

        t = time.time()
        features = representation.transform(self.frames[frame_idx])

        verbosity_wrap('This took {} seconds/frame.'.format(round((time.time()-t)/len(frame_idx),8)))

        if(pretend==False):
            verbosity_wrap("Estimating time to compute kernel...")
            self.estimate_time(N=max(0,20-len(self.est_frames)))
            est = int(np.poly1d(np.polyfit(self.est_frames,self.est_times,deg=2))(len(frame_idx)))+1
        else:
            est=0

        verbosity_wrap("<br/>Next we find the kernel for our training model.<br/>(This step will take approximately {} minutes and {} seconds.)".format(int(est/60), int(est%60)))
        time.sleep(0.5)
        t = time.time()
        kernel = compute_kernel(representation, self.zeta, features, kernel_type=self.sliders['kernel_type'].value, time_estimate=est)
        self.est_frames.append(len(frame_idx))
        self.est_times.append(time.time()-t)

        delta = np.std(training_properties) / np.mean(kernel.diagonal())
        kernel[np.diag_indices_from(kernel)] += self.Lambda**2 / delta **2 + jitter

        verbosity_wrap("<br/>We will adjust the diagonals of our kernel by {} so that it is properly scaled.".format(round(self.Lambda**2 / delta **2 + jitter,8)))
        verbosity_wrap("<br/>Now we can take this kernel to compute the weights of our KRR.")

        weights = np.linalg.solve(kernel,training_properties)
        model = KRR(self.zeta, weights, representation, features, kernel_type=self.sliders['kernel_type'].value)
        if(pretend==False):
            self.krr[self.sliders['property_to_ml'].value], k = model, kernel
            self.trained[self.sliders['property_to_ml'].value] = True
    def predict_test_set(self):
        self.plot_prediction_func(frame_idx = self.test_idx)
    def predict_new_set(self, filename = './reference_data/small_molecules-1000.xyz', num_frames=''):
        frames=np.array(read(filename, ":{}".format(num_frames)))
        properties = extract_property(frames, self.sliders['property_to_ml'].value)
        self.plot_prediction_func(y_known=properties, frames=frames)
    def plot_prediction_func(self, y_known=None, frames=None, frame_idx=[], pretend=False):
        if(len(frame_idx)>0):
            frames=self.frames[frame_idx]
            y_known = np.concatenate(self.properties[self.sliders['property_to_ml'].value][frame_idx])[:,0] if self.sliders['kernel_type'].value=='atomic' else self.properties[self.sliders['property_to_ml'].value][frame_idx]
        if(pretend==False):
            verbosity_wrap = lambda s: self.verbosity_wrap(s)
        else:
            verbosity_wrap = lambda s: None
        if(self.trained[self.sliders['property_to_ml'].value]==False):
            verbosity_wrap("Model has not yet been trained, training now..."   )
            self.train_krr_model()
        testing_properties = y_known#np.concatenate(self.properties[self.sliders['property_to_ml'].value][frame_idx])[:,0] if self.sliders['kernel_type'].value=='atomic' else self.properties[self.sliders['property_to_ml'].value][frame_idx]

        if(pretend==False):
            verbosity_wrap("Estimating time to compute prediction...")
            self.estimate_time(x=self.pred_frames, y=self.pred_times, f=self.plot_prediction_func, N=max(0,20-len(self.pred_frames)), ref=len(frames))
            est = int(np.poly1d(np.polyfit(self.pred_frames,self.pred_times,deg=3))(len(frames)))+1
        else:
            est=0

        t = time.time()
        verbosity_wrap("Predicting the properties of our data set will take approximately {} minutes and {} seconds.".format(int(est/60), int(est%60)))
        y_pred = self.krr[self.sliders['property_to_ml'].value].predict(frames)
        self.pred_frames.append(len(frame_idx))
        self.pred_times.append(time.time()-t)
        verbosity_wrap(markdown_table_from_dict(get_score(y_pred, testing_properties), headers=["Statistic","Value" ]))

        if(not pretend):
            plt.scatter(y_pred, testing_properties, s=3)
            plt.axis('scaled')
            plt.xlabel(self.sliders['property_to_ml'].value)
            plt.ylabel('Predicted '+self.sliders['property_to_ml'].value)
            plt.gca().set_aspect('equal')
            plt.show()
    def output_params(self):
        self.verbosity_wrap('Our input file is {}, of which we are using {} frames.'.format(self.input_file, self.sliders['number_of_frames'].value))
        self.verbosity_wrap("<br/>Our hyperparameters are {}".format(markdown_table_from_dict({hyper_vals[k]['name']:[self.hyperparameters[k]] for k in self.hyperparameters if 'fixed' not in hyper_vals[k]}, headers =["Parameter", "Value"])))
    def set(self, value_name, value):
        assert value_name in self.sliders
        self.sliders[value_name].value=value
    def estimate_time(self, x=None, y=None, f=None, N=20, ref =None):
        if(N<=0):
            return
        if(x==None or y==None or f==None or ref==None):
            x=self.est_frames
            y=self.est_times
            f=self.train_krr_model_func
            ref = int(0.25*self.sliders['number_of_frames'].value*self.sliders['training_percentage'].value)

        for nf in tqdm(np.random.randint(2,ref, size=min(N, ref))):
            f(frame_idx=np.random.randint(len(self.frames),size=nf), pretend=True)
