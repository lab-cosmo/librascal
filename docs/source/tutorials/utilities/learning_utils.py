import os
import inspect
import time

from IPython.display import Markdown, display, clear_output, Image, HTML
import ipywidgets as widgets
import numpy as np
import ase
from ase.io import read
from matplotlib import pyplot as plt
from tqdm import tqdm_notebook as tqdm
import pandas as pd
import json

# This is for backwards compatibility -- old versions of RASCAL will follow
# the latter schema, newer versions the former.
try:
    from rascal.representations import SphericalInvariants as SOAP
except:
    from rascal.representations import SOAP

from general_utils import (markdown_table_from_dict, \
                           readme_button, _button_template_,
                           compute_kernel, mask_body_order)

# These are some imports/helpers for all of the tutorials
style = json.load(open('./utilities/widget_style.json'))
hyper_dict = json.load(open("./utilities/hyperparameter_presets.json"))
hyper_vals = json.load(open("./utilities/hyperparameter_bounds.json"))

# Properties preset for the learning tutorial
known_properties = dict(CS="Atom",
                        dft_formation_energy_per_atom_in_eV="Structure")
ignore = ['Natoms', 'numbers', 'Name',
    'positions', 'cutoff', 'nneightol', "NAME"]

def split_dataset(frames, train_fraction, seed=10):
    """ Function to split a dataset into testing and training portions.

        Author: M. Veit or F. Musil

        Keywords:
            frames - ASE-Atoms like
            train_fraction - float in [0.0,1.0], percentage of frames to train

        Returns:
            indices of training and testing frames
    """

    N = len(frames)
    ids = np.arange(N)
    np.random.seed(seed)
    np.random.shuffle(ids)
    Ntrain = int(N*train_fraction)
    train = ids[:Ntrain]
    test = ids[Ntrain:]
    return np.array(train), np.array(test)


def get_r2(y_pred, y_true):
    """ Function to calculate R^2
        Author: M. Veit or F. Musil

        Keywords:
            y_pred - list of predicted y values
            y_true - list of known y values

        Returns:
            float
    """
    numerator = ((y_true - y_pred) ** 2).sum(axis=0, dtype=np.float64)
    denominator = ((y_true - np.average(y_true,
                                        axis=0)) ** 2).sum(axis=0,
                                                            dtype=np.float64)
    output_scores = 1 - numerator / denominator
    return np.mean(output_scores)


def get_score(y_pred, y_true):
    """ Function to calculate accuracy statistics

        Author: M. Veit or F. Musil

        Keywords:
            y_pred - list of predicted y values
            y_true - list of known y values

        Returns:
            dictionary containing SUP (supremum), MAE (mean squared error),
                                  RMSD (root mean squared displacement),
                                  and R^2
    """
    return {
                "SUP": [np.amax(np.abs((y_pred-y_true)))],
                "MAE": [np.mean(np.abs(y_pred-y_true))],
                "RMSD": [np.sqrt(np.mean((y_pred-y_true)**2))],
                r"$R^2$": [get_r2(y_pred, y_true)]
                }


class KRR(object):
    """ Class for Kernel Ridge Regression

    Author: M. Veit or F. Musil

    Keywords:
        zeta - integer
        weights - weights of the given by the kernel and the targets given by
                  the training set
        representation - feature set given by soap vectors
        kernel_type - "atom" or "Structure"
        X - training dataset

    Functions:
        predict - given the KRR, predicts the properties for a set of ASE frames
    """

    def __init__(self, weights, features, kernel_type, calculator=None, **kwargs):
        self.weights = weights
        self.hypers = dict(**kwargs)
        self.calculator = calculator if calculator!=None else SOAP(**mask_body_order(kwargs))
        self.X = features
        self.kernel_type = kernel_type

    def predict(self, frames):
        features = self.calculator.transform(frames)
        kernel = compute_kernel(calculator=self.calculator, features1=self.X, \
                                features2=features, kernel_type=self.kernel_type, \
                                **self.hypers)
        return np.dot(self.weights, kernel)


def extract_property(frames, property='energy'):
    """ Function to read properties from ASE frames object

    Author: M. Veit or F. Musil

    Keywords:
        frames - ASE Atoms-like object
        property - property to extract

    Returns:
        list of the property from the given frames

    Raises error when property is not found in the ASE object
    """
    if(property in frames[0].info):
        return np.array([cc.info[property] for cc in frames])
    elif(property in frames[0].arrays):
        return np.array([cc.arrays[property] for cc in frames])
    else:
        print(frames[0].info)
        raise KeyError(
            "{} is not a property in the given frames".format(property))


class learning_tutorial(object):
    """ Class to wrap the learning tutorial

    Author: R. Cersonsky

    Constructor Keywords:
        input_file          - ASE Atoms-like object containing molecules or
                              crystals with some property to learn
        training_percentage - float, percentage of frames to include in the
                              training set
        interactive         - boolean, whether to load the tutorial interactive-
                              ly in the jupyter notebook
        verbose             - boolean, whether or not to narrate the actions of
                              the tutorial
        hyperparameters     - dictionary, hyperparameters for SOAP calculation
        number_of_frames    - int, how many frames to include in the dataset
        property            - string, which property to train the dataset on

    """

    def __init__(self, input_file='./data/learning/small_molecules-1000.xyz',
                 training_percentage=0.8, interactive=False, verbose=True,
                 hyperparameters=dict(**hyper_dict['Power Spectrum']),
                 number_of_frames=None, property=None):

        # presets given by M. Veit and F. Musil
        self.zeta = 2
        self.Lambda = 5e-3

        self.hyperparameters = hyperparameters
        self.verbose = verbose

        # The verbosity wrap serves to limit the amount of text output in the
        # jupyter notebook
        self.verbosity_wrap = lambda s: (None if not verbose
                                              else display(Markdown(s)))

        # Looks in the examples folder for files suitable for the learning
        # tutorial
        file_options = ['./data/learning/{}'.format(f) for f in
                            os.listdir('./data/learning/') if f.endswith('xyz')]

        # If the supplied file is not in the examples folder, includes it in the
        # list of suitable files
        if(input_file != None):
            file_options.insert(0, input_file)
            if(input_file in file_options[1:]):
                file_options.pop(file_options[1:].index(input_file)+1)

        # Sets the file to the default if not supplied
        self.input_file = file_options[0] if input_file == None else input_file

        # Reads input file
        self.frames = np.array(read(self.input_file, ":"))
        self.train_idx, self.test_idx = None, None

        # This next section is for the jupyter interface. Even if the user has
        # specified zero interactivity, the sliders serve as containers for the
        # hyperparameters to change when set
        self.sliders = {val:
                        widgets.FloatSlider(
                            value=self.hyperparameters.get(val,
                                                hyper_vals[val]['options'][0]),
                            min=hyper_vals[val]['options'][0],
                            max=hyper_vals[val]['options'][1],
                            description=hyper_vals[val]['name'],
                            continuous_update=True,
                            step=(hyper_vals[val]['options'][1] -
                                    hyper_vals[val]['options'][0])/20.,
                            style=style,
                            )
                        if isinstance(hyper_vals[val]['options'][0], float) else
                        widgets.IntSlider(
                            value=self.hyperparameters.get(val,
                                                hyper_vals[val]['options'][0]),
                            min=hyper_vals[val]['options'][0],
                            max=hyper_vals[val]['options'][1],
                            description=hyper_vals[val]['name'],
                            continuous_update=True,
                            step=1,
                            style=style,
                            )
                        if isinstance(hyper_vals[val]['options'][0], int)
                            and hyper_vals[val]['options'][0] != True
                        else widgets.Dropdown(
                                options=hyper_vals[val]['options'],
                                style=style,
                                value=self.hyperparameters.get(val,
                                                hyper_vals[val]['options'][0]),
                                description=hyper_vals[val]['name'],
                                )
                        for val in hyper_vals
                        if 'fixed' not in hyper_vals[val]}
        self.properties = {prop:
                            extract_property(self.frames, prop)
                            for prop in sorted([*list(self.frames[0].info.keys()),
                                                *list(self.frames[0].arrays.keys())],
                                                key=lambda p: -int(p in known_properties))
                            if prop not in ignore}

        self.sliders['number_of_frames'] = widgets.IntSlider(
                                                value=int(len(self.frames)*0.2)
                                                      if number_of_frames == None
                                                      else number_of_frames,
                                                min=1,
                                                max=len(self.frames),
                                                description="Number of Frames",
                                                step=1,
                                                style=style)
        self.sliders['property_to_ml'] = widgets.Dropdown(
                                                value=list(
                                                    self.properties.keys())[0]
                                                      if property == None
                                                      else property,
                                                options=list(
                                                    self.properties.keys()),
                                                description="Property to ML",
                                                style=style)
        if(self.sliders['property_to_ml'].value not in known_properties):
            self.sliders['kernel_type'] = widgets.Dropdown(
                                                value='Atom',
                                                options=['Atom', 'Structure'],
                                                description="Kernel Type",
                                                style=style)
        else:
            self.sliders['kernel_type'] = widgets.Dropdown(
                                                value=known_properties[self.sliders['property_to_ml'].value],
                                                options=[
                                                    known_properties[self.sliders['property_to_ml'].value]],
                                                description="Kernel Type",
                                                style=style)
        self.sliders['training_percentage'] = widgets.FloatSlider(
                                                    value=training_percentage,
                                                    min=0,
                                                    max=1,
                                                    description="Training Percentage",
                                                    continuous_update=True,
                                                    step=0.05,
                                                    style=style,
                                                    )

        # Show and enable the widgets if interactive==True
        if(interactive):
            self.input_button = widgets.Dropdown(options=[*file_options],  # , "Other"],
                                         style=style,\
                                         value=self.input_file,
                                         description="Input File: ",
                                         )
            self.input_button.observe(self.get_input, names='value')
            display(self.input_button)

            self.preset_button = _button_template_(list(hyper_dict.keys()),
                                                "SOAP Presets: ", self.preset_func)

            slider_order = ['property_to_ml', 'kernel_type',
                            'number_of_frames', 'training_percentage',
                            *list(hyper_vals.keys())]
            for s in slider_order:
                if(s in self.sliders):
                    self.sliders[s].observe(lambda change:
                                            self.change_func(change),
                                            names='value')
                    display(self.sliders[s])

        # Set up one KRR for each possible property for caching reasons
        self.krr = {prop: None for prop in self.properties}
        self.trained = {prop: False for prop in self.properties}

        # These serve as placeholders for time estimates until any calculation
        # is run
        self.est_frames, self.est_times = [], []
        self.pred_frames, self.pred_times = [], []

    def reset_ML(self, inp_change=False):
        """ Function to reset the properties, number of frames, and kernel type
            if any values affecting the machine learning have changed

            Author: R. Cersonsky
        """
        if(inp_change):
            self.properties = {prop:
                extract_property(self.frames, prop)
                for prop in sorted([*list(self.frames[0].info.keys()),
                                    *list(self.frames[0].arrays.keys())],
                                    key=lambda p: -int(p in known_properties))
                if prop not in ignore}

            self.sliders['number_of_frames'].max = len(self.frames)
            self.sliders['number_of_frames'].value = int(len(self.frames)*0.2)

            self.sliders['property_to_ml'].options = list(
                self.properties.keys())
            self.sliders['property_to_ml'].value = list(
                self.properties.keys())[0]

            if(self.sliders['property_to_ml'].value in known_properties):
                self.sliders['kernel_type'].options = [
                    known_properties[self.sliders['property_to_ml'].value]]
                self.sliders['kernel_type'].value = self.sliders['kernel_type'].options[0]
            else:
                self.sliders['kernel_type'].options = ['Atom', 'Structure']
                self.sliders['kernel_type'].value = 'Atom'

        self.krr = {prop: None for prop in self.properties}
        self.trained = {prop: False for prop in self.properties}

    def change_func(self, change):
        """ Function for button events

            Author: R. Cersonsky
        """
        change['owner'].value = change['new']
        for s in self.hyperparameters:
            if(s in self.sliders):
                if(self.hyperparameters[s] != self.sliders[s].value):
                    self.est_frames, self.est_times = [], []
                    self.hyperparameters[s] = self.sliders[s].value

        if(change['owner'] == self.sliders['property_to_ml']):
            if(self.sliders['property_to_ml'].value in known_properties):
                self.sliders['kernel_type'].options = [
                    known_properties[self.sliders['property_to_ml'].value]]
                self.sliders['kernel_type'].value = self.sliders['kernel_type'].options[0]
            else:
                self.sliders['kernel_type'].options = ['Atom', 'Structure']
                self.sliders['kernel_type'].value = 'Atom'

        self.reset_ML()

    def preset_func(self, a):
        """ Function to change hyperparameters based upon presets

            Author: R. Cersonsky
        """
        for s in self.sliders:
            if(s in self.hyperparameters):
                self.hyperparameters[s] = hyper_dict[self.preset_button.value][s]
                self.sliders[s].value = self.hyperparameters[s]
        self.reset_ML()

    def get_input(self, a):
        """ Function to change input file

            Author: R. Cersonsky
        """
        inp = self.input_button.value
        self.input_file = inp
        self.frames = np.array(read(self.input_file, ":"))
        self.reset_ML(True)

    def train_krr_model_func(self, frame_idx, jitter=1e-8, pretend=False):
        """ Function to trains the model on the given SOAP representation and
            properties. It is used for both the time estimation and the eventual
            training. The flag pretend is set to True when it is being used for
            estimation, false otherwise.

            Author: R. Cersonsky

            Keywords:
                frame_idx - list of indices of the frames to train on
                jitter    - float, anticipated error to adjust kernel diagonals
                pretend   - boolean, whether to not show the results

        """

        if(pretend == False):
            self.output_params()
            verbosity_wrap = lambda s: self.verbosity_wrap(s)
        else:
            verbosity_wrap = lambda s: None

        verbosity_wrap("<br/>We will now train a model on {}.".format(
                                        self.sliders['property_to_ml'].value))

        representation = SOAP(**mask_body_order(self.hyperparameters))

        if(known_properties[self.sliders['property_to_ml'].value] == 'Atom'):
            props = np.concatenate(self.properties[self.sliders['property_to_ml'].value][frame_idx])[:,0]
        else:
            props = self.properties[self.sliders['property_to_ml'].value][frame_idx]

        training_properties=props

        verbosity_wrap("First, I am going to separate my dataset:")
        verbosity_wrap("<br/>Now we will compute the SOAP representation of \
                        our training frames.")

        t=time.time()
        features=representation.transform(self.frames[frame_idx])

        verbosity_wrap('This took {} seconds/frame.'.format(
                                    round((time.time()-t)/len(frame_idx), 8)))

        if(pretend == False):
            self.estimate_time(N=max(0, 20-len(self.est_frames)),
                               p="Estimating time to compute kernel...")
            est=int(np.poly1d(np.polyfit(self.est_frames,
                                           self.est_times,
                                           deg=2))(len(frame_idx)))+1
        else:
            est=0

        verbosity_wrap("<br/>Next we find the kernel for our training model.\
                        <br/>(This step will take approximately {} minutes and \
                        {} seconds.)".format(int(est/60), int(est % 60)))

        time.sleep(0.5)
        t=time.time()
        kernel = compute_kernel(calculator=representation, \
                                features1=features, \
                                kernel_type=self.sliders['kernel_type'].value, \
                                **self.hyperparameters)

        self.est_frames.append(len(frame_idx))
        self.est_times.append(time.time()-t)

        delta=np.std(training_properties) / np.mean(kernel.diagonal())
        kernel[np.diag_indices_from(
            kernel)] += self.Lambda**2 / delta ** 2 + jitter

        verbosity_wrap("<br/>We will adjust the our kernel with the tolerance \
                        matrix $\\Lambda = ({})I$.".format(\
                        round(self.Lambda**2 / delta ** 2 + jitter, 8)))
        verbosity_wrap("<br/>Now we can take this kernel to compute the weights \
                        of our KRR.")

        weights=np.linalg.solve(kernel, training_properties)
        model=KRR(weights, features, \
                  kernel_type=self.sliders['kernel_type'].value,
                  **self.hyperparameters
                  )
        if(pretend == False):
            self.krr[self.sliders['property_to_ml'].value], k=model, kernel
            self.trained[self.sliders['property_to_ml'].value]=True

    def train_krr_model(self, jitter= 1e-8):
        """ Wrapper to the train_krr_model_func function for use in the tutorial

            Author: R. Cersonsky

            Keywords:
                jitter    - float, anticipated error to adjust kernel diagonals

        """

        training_dict={"Training Set": [int(self.sliders['number_of_frames'].value*(self.sliders['training_percentage'].value)),\
                                          round(100*(self.sliders['training_percentage'].value))],
                         "Testing Set": [int(self.sliders['number_of_frames'].value*(1-self.sliders['training_percentage'].value)),\
                                          round(100*(1-self.sliders['training_percentage'].value))]}
        headers=["Partition", "Number of Frames", "Percentage"]
        self.verbosity_wrap(markdown_table_from_dict(training_dict, headers))

        self.train_idx, self.test_idx=split_dataset(self.frames[: self.sliders['number_of_frames'].value],
                                                      self.sliders['training_percentage'].value)
        self.train_krr_model_func(self.train_idx, jitter= jitter)

    def plot_prediction_func(self, y_known= None, frames = None,
                             frame_idx= [], pretend = False):
        """ Function to predict properties from the given SOAP representation and
            KRR model. It is used for both the time estimation and the tutorial
            prediction. The flag pretend is set to True when it is being used for
            estimation, false otherwise.

            Author: R. Cersonsky

            Keywords:
                y_known   - list, target properties to train on
                frames    - frames to predict properties of
                frame_idx - list of indices of the frames to train on
                pretend   - boolean, whether to not show the results
        """
        if(len(frame_idx) > 0):
            frames=self.frames[frame_idx]
            y_known = np.concatenate(self.properties[self.sliders['property_to_ml'].value][frame_idx])[: , 0] if self.sliders['kernel_type'].value == 'Atom' else self.properties[self.sliders['property_to_ml'].value][frame_idx]
        if(pretend == False):
            verbosity_wrap=lambda s: self.verbosity_wrap(s)
        else:
            verbosity_wrap = lambda s: None
        if(self.trained[self.sliders['property_to_ml'].value] == False):
            verbosity_wrap("Model has not yet been trained, training now...")
            self.train_krr_model()
        # np.concatenate(self.properties[self.sliders['property_to_ml'].value][frame_idx])[:,0] if self.sliders['kernel_type'].value=='Atom' else self.properties[self.sliders['property_to_ml'].value][frame_idx]
        testing_properties=y_known

        if(pretend == False):
            self.estimate_time(x=self.pred_frames, y=self.pred_times, f=self.plot_prediction_func, N=max(
                0, 20-len(self.pred_frames)), ref=len(frames), p="Estimating time to compute prediction...")
            est=int(np.poly1d(np.polyfit(self.pred_frames,
                    self.pred_times, deg=3))(len(frames)))+1
        else:
            est=0

        t=time.time()
        verbosity_wrap("Predicting the properties of our data set will take approximately {} minutes and {} seconds.".format(
            int(est/60), int(est % 60)))
        y_pred=self.krr[self.sliders['property_to_ml'].value].predict(frames)
        self.pred_frames.append(len(frame_idx))
        self.pred_times.append(time.time()-t)
        verbosity_wrap(markdown_table_from_dict(
            get_score(y_pred, testing_properties), headers = ["Statistic", "Value"]))

        if(not pretend):
            plt.scatter(y_pred, testing_properties, s=3)
            plt.axis('scaled')
            plt.xlabel(self.sliders['property_to_ml'].value)
            plt.ylabel('Predicted '+self.sliders['property_to_ml'].value)
            plt.gca().set_aspect('equal')
            plt.show()

    def predict_test_set(self):
        """ Wrapper for plot_prediction_func that is used within the tutorial

            Author: R. Cersonsky
        """
        self.plot_prediction_func(frame_idx=self.test_idx)

    def predict_new_set(self, filename='./data/small_molecules-1000.xyz',
                        num_frames=''):
        """ Wrapper predicts a new set of data beyond the training set, loaded
            from a different file

            Author: R. Cersonsky

            Keywords:
                filename   - string corresponding to an ASE-like file
                num_frames - number of frames to predict
        """
        frames=np.array(read(filename, ":{}".format(num_frames)))
        properties=extract_property(
            frames, self.sliders['property_to_ml'].value)
        self.plot_prediction_func(y_known=properties, frames=frames)

    def output_params(self):
        """ Function to output the hyperparameters that are currently set

            Author: R. Cersonsky
        """
        self.verbosity_wrap('Our input file is {}, of which we are using {} \
                            frames.'.format(self.input_file, \
                                            self.sliders['number_of_frames'].value))

        hdict={hyper_vals[k]['name']: [self.hyperparameters[k]]
                                                for k in self.hyperparameters
                                                if 'fixed' not in hyper_vals[k]}
        headers=["Parameter", "Value"]

        self.verbosity_wrap("<br/>Our hyperparameters are {}".format(
                                      markdown_table_from_dict(hdict, headers)))

    def set(self, value_name, value):
        """ Public Setter function for tutorial interface

            Author: R. Cersonsky

            Keywords:
            value_name - string, must be in self.sliders
            value      - value to set
        """
        assert value_name in self.sliders
        self.sliders[value_name].value=value

    def estimate_time(self, x=None, y=None, f=None, N=20, ref=None, p=None):
        """ Helper function for estimating time

            Author: R. Cersonsky

            Keywords:
            x   - list of numbers of frames computed
            y   - list, size=size(x), of times to computed frames in x
            f   - function to estimate the run time of
            N   - number of trials to run
            ref - int, upper limit of trials to run.
                  Defaults to (number_of_frames in the ultimate computation) / 4
            p   - print statement
        """
        if(N <= 0):
            return
        if(x == None or y == None or f == None or ref == None):
            x=self.est_frames
            y=self.est_times
            f=self.train_krr_model_func
            ref=int(0.25*self.sliders['number_of_frames'].value *
                    self.sliders['training_percentage'].value)

        if(p != None and min(N, ref) > 0):
            self.verbosity_wrap(p)
        for nf in tqdm(np.random.randint(2, ref, size = min(N, ref))):
            f(frame_idx=np.random.randint(len(self.frames), size = nf), pretend=True)
