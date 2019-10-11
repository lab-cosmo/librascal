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


def format_positions(positions):
    """ Function to format positions for best viewing angle. Does a simple PCA
        to find the best direction to view the molecule and rotates position
        along these axes

    Author: R. Cersonsky

    Keywords:
        positions - list of 3-vectors

    Returns:
        positions - list of 3-vectors
    """
    from sklearn.decomposition import PCA
    pca = PCA(n_components=3)
    return pca.fit_transform(positions)


class SOAP_tutorial(object):
    """ Class to wrap the SOAP vectors tutorial

    Author: R. Cersonsky

    Constructor Keywords:
        molecule_file       - filename of ASE Atoms-like object
                              containing molecule
        interactive         - boolean, whether to load the tutorial interactive-
                              ly in the jupyter notebook
        verbose             - boolean, whether or not to narrate the actions of
                              the tutorial
        hyperparameters     - dictionary, hyperparameters for SOAP calculation
                              or string corresponding to one of the SOAP presets

    Functions:
        predict - given the KRR, predicts the properties for a set of ASE frames
    """

    def __init__(self,
                 molecule_file=None,
                 interactive=False,
                 verbose=True,
                 hyperparameters=dict(**hyper_dict['Power Spectrum']),
                 ):

        # Check if hyperparameters are of type string or dictionary
        if(isinstance(hyperparameters, str)):
            self.hyperparameters = hyper_dict[hyperparameters]
        else:
            self.hyperparameters = {h: hyperparameters[h]
                                    if h in hyperparameters
                                    else hyper_dict['Power Spectrum'][h]
                                    for h in hyper_dict['Power Spectrum']}
        self.verbose = verbose
        self.verbosity_wrap = lambda s: (
            None if not verbose else display(Markdown(s)))

        # Look for and load available molecules
        file_options = ['./data/molecules/{}'.format(f)
                        for f in os.listdir('./data/molecules/')
                        if f.endswith('xyz')]

        if(molecule_file != None and molecule_file not in file_options):
            file_options.append(molecule_file)

        self.frames, self.positions, self.symbols, self.inds, self.labels = {}, {}, {}, {}, {}
        self.import_molecules(file_options)

        self.vectors = None
        self.interactive = interactive

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
                        if isinstance(hyper_vals[val]['options'][0], int) and
                        hyper_vals[val]['options'][0] != True else
                        widgets.Dropdown(options=hyper_vals[val]['options'],
                                         style=style,
                                         value=self.hyperparameters.get(val,
                                                                        hyper_vals[val]['options'][0]),
                                         description=hyper_vals[val]['name'],
                                         )
                        for val in hyper_vals if 'fixed' not in hyper_vals[val]}

        # Which elements to center smooth gaussians on. Defaults to all elements
        # present except for hydrogen.
        constituents = list(set([s for k in self.symbols for s in self.symbols[k]]))
        self.sliders['center_select']=widgets.SelectMultiple(
                                    options=constituents,
                                    value=[s for s in constituents if s != "H"],
                                    description="Where to Place SOAP Centers",
                                    style=style)

        # Whether to compute the global kernel (average of the environments) or
        # atomic kernel.
        self.sliders['average'] = widgets.Dropdown(
            options=["Environment-Centered", "Average"],
            value="Environment-Centered",
            description="Type of SOAP Vectors",
            style=style)

        # Show and enable the widgets if interactive==True
        if(interactive):
            self.preset_button = _button_template_(list(hyper_dict.keys()),
                                                   "SOAP Presets: ",
                                                   self.preset_func)
            if(isinstance(hyperparameters, str)):
                self.preset_button.value = hyperparameters

            slider_order = ['center_select',
                            'average',
                            *list(hyper_vals.keys())]

            for s in slider_order:
                if(s in self.sliders):
                    display(self.sliders[s])
                    self.sliders[s].observe(
                        lambda change: self.change_func(change), names='value')

        # These buttons allow users to choose up to two molecules to compute
        # SOAP vectors and compare.
        self.molecule_buttons = [widgets.Dropdown(options=[*self.frames.keys()],  # , "Other"],
                                                  style=style,\
                                                  value=list(
                                                      self.frames.keys())[0],
                                                  description="\tMolecule: "
                                                  ) for i in range(2)]

    def import_molecules(self, file_list):
        """ Function to import ASE-like files and include them in the dictionaries
            of available molecules.

            Author: R. Cersonsky

            Keywords:
                 file_list - list of strings corresponding to ASE-like files
        """
        for file in file_list:
            for frame in read(file, ":"):
                key = str(frame.symbols)

                if(key in self.frames):
                    print("Overwriting previous entry for {}".format(key))

                self.frames[key] = frame
                self.positions[key] = format_positions(frame.get_positions())
                self.inds[key] = sorted(range(len(self.positions[key])),
                                        key=lambda j: self.positions[key][j][0])
                self.positions[key] = self.positions[key][self.inds[key]]
                self.symbols[key] = frame.symbols[self.inds[key]]
                super_scripts = ['' if list(self.symbols[key]).count(s) == 1
                                 else '^{{({})}}'.format(list(self.symbols[key])[:i].count(s)+1)
                                 for i, s in enumerate(list(self.symbols[key]))]
                self.labels[key] = np.array([r'${}{}$'.format(self.symbols[key][i], s)
                                             for i, s in enumerate(super_scripts)])

    def change_func(self, change):
        """ Function to import ASE-like files and include them in the dictionaries
            of available molecules.

            Author: R. Cersonsky

            Keywords:
                 file_list - list of strings corresponding to ASE-like files
        """
        change['owner'].value = change['new']
        for s in self.hyperparameters:
            if(s in self.sliders):
                if(self.hyperparameters[s] != self.sliders[s].value):
                    self.hyperparameters[s] = self.sliders[s].value
                    self.vectors = None

    def preset_func(self, a):
        """ Function to change hyperparameters based upon presets

            Author: R. Cersonsky
        """
        for s in self.sliders:
            if(s in self.hyperparameters):
                self.hyperparameters[s] = hyper_dict[self.preset_button.value][s]
                self.sliders[s].value = self.hyperparameters[s]

    def get_soap_vectors(self, key=None, average=False):
        """ Get the SOAP vectors for the given molecule, denoted by key

            Author: R. Cersonsky

            Keywords:
                 key     - string corresponding to molecule built by \
                           import_molecules
                 average - whether to compute the average SOAP vector for the
                           molecule
        """
        if(key == None):
            key = self.input_file.replace('./data/molecules/', '')[:-4]

        representation = SOAP(**mask_body_order(self.hyperparameters))
        all_frames = [self.frames[k] for k in sorted(self.frames.keys())]
        features = representation.transform(
            all_frames).get_dense_feature_matrix(representation)

        lengths = [len(self.inds[k]) for k in sorted(self.frames.keys())]

        self.vectors = {k: features[sum(lengths[0:i]):sum(lengths[0:i])+lengths[i]]
                        for i, k in enumerate(sorted(self.frames.keys()))}
        vectors = self.vectors[key]

        # if(average):
        #     return np.mean(vectors, axis=0)
        return vectors[self.inds[key]]

    def make_figure(self, key, img_name):
        """ Function to generate an image of the given molecule and cache it
            under img_name

            Author: R. Cersonsky

            Keywords:
                 key      - string corresponding to molecule built by
                            import_molecules
                 img_name - filename under which to cache the image
        """
        from ase.data import covalent_radii, atomic_numbers
        from matplotlib import pyplot, patches
        from ase.data.colors import jmol_colors
        style = {
            "horizontalalignment": "center",
            "verticalalignment": "center",
            "fontsize": 12
        }
        fig, ax = pyplot.subplots(1)

        for i, p in enumerate(self.positions[key]):
            an = atomic_numbers[self.symbols[key][i]]
            ax.text(*p[:2],
                    s=self.labels[key][i],
                    **style,
                    zorder=p[-1]+0.01
                    )
            ax.add_artist(patches.Circle(p[:2],
                                         covalent_radii[an],
                                         facecolor=jmol_colors[an],
                                         edgecolor='k',
                                         alpha=0.95,
                                         zorder=p[-1]))
        ax.axis('off')
        bounds = [*np.min(self.positions[key], axis=0), *
                  np.max(self.positions[key], axis=0)]
        ax.set_xlim([bounds[0]-1.5, bounds[3]+1.5])
        ax.set_ylim([bounds[1]-1.5, bounds[5]+1.5])
        ax.set_aspect('equal')
        fig.savefig(img_name)
        fig.clf()

    def show_kernel(self, key1, key2=None, show=True, average=False):
        from rascal.models import Kernel
        from ase.data import chemical_symbols
        """ Function to show pandas dataframe for the comparison kernel

            Author: R. Cersonsky

            Keywords:
                key1    - string corresponding to molecule built by
                          import_molecules
                key2    - string corresponding to molecule built by
                          import_molecules, defaults to key1
                show    - boolean, whether to show or return the dataframe
                average - boolean, whether to compute the average SOAP vector
                          for the molecule
        """
        key2 = key2 if key2 != None else key1
        average = self.sliders['average'].value == "Average"

        soap=SOAP(**mask_body_order(self.hyperparameters))
        features1=soap.transform(self.frames[key1])
        features2=soap.transform(self.frames[key2])

        kernel = Kernel(soap,
                        target_type="Structure" if average else "Atom",
                        zeta = 2,
                        **self.hyperparameters)
        data = kernel(features2, features1)

        if (average):

            vec1 = np.mean(features1.get_dense_feature_matrix(soap),axis=0)
            vec2 = np.mean(features2.get_dense_feature_matrix(soap), axis=0)
            if(vec1.shape==vec2.shape):
                data = np.dot(vec1, vec2)/(np.linalg.norm(vec1)*np.linalg.norm(vec2))
            df=pd.DataFrame(data=[['%1.3f' %(data)]],
                            columns =[key2],
                            index =[key1]
                        )
        else:
            df=pd.DataFrame(data=[['%1.3f' %(data[i][j])
                                    for i,v in enumerate(features2._frames.numbers) if v!=1]
                                    for j,w in enumerate(features1._frames.numbers) if w!=1],
                            columns =[chemical_symbols[v] for i,v in enumerate(features2._frames.numbers) if v!=1],
                            index =[chemical_symbols[v] for i,v in enumerate(features1._frames.numbers) if v!=1]
                            )
        #
        df.name = ("Self-" if key1 == key2 else '') + \
            "Similarity " + ("Kernel" if not average else "")
        return df

    def show_molecule(self, key, show=True):
        """ Function to either generate or show the image for the given molecule

            Author: R. Cersonsky

            Keywords:
                key  - string corresponding to molecule built by
                       import_molecules
                show - boolean, whether to show or return the dataframe
        """
        import os
        img_name = './images/cache/{}.png'.format(key)
        if(not os.path.exists(img_name)):
            self.make_figure(key, img_name)
        if(show):
            display(Image(img_name))
        return img_name

    def compute(self, a=None):
        """ Function to compute the SOAP vectors for the set hyperparameters

            Author: R. Cersonsky
        """
        def wrap_compute(a):
            self.get_soap_vectors(
                key=self.molecule_buttons[0].value, average=self.sliders['average'].value == "Average")
            print("SOAP Vectors for {}: \n\t".format(
                self.molecule_buttons[0].value), self.vectors[self.molecule_buttons[0].value])

        compute_button = widgets.Button(description="Compute",
                                        style=style,
                                        value=False,
                                        )
        compute_button.on_click(wrap_compute)
        clear_button = widgets.Button(description="Clear Output",
                                      style=style,
                                      value=False,
                                      )
        clear_button.on_click(lambda a: self.compute(a))
        display(self.molecule_buttons[0])
        display(widgets.widgets.HBox((compute_button, clear_button)))

    def compare(self, a=None):
        """ Function to compare the SOAP vectors of two molecules
            for the set hyperparameters

            Author: R. Cersonsky
        """
        clear_output()

        def wrap_compare(a):
            df = self.show_kernel(
                key1=self.molecule_buttons[0].value, key2=self.molecule_buttons[1].value, show=False)
            m1 = self.show_molecule(self.molecule_buttons[0].value, show=False)
            m2 = self.show_molecule(self.molecule_buttons[1].value, show=False)
            display(HTML("<table style='align: center'><tr><th style=\"text-align:center\">{}</th><th style=\"text-align:center\">{}</th><th style=\"text-align:center\">{}</th><tr><td><img src={}></td><td>{}</td><td><img src={}></td></tr></table>".format(
                self.molecule_buttons[0].value, df.name, self.molecule_buttons[1].value, m1, df.to_html(), m2)))
        compare_button = widgets.Button(description="Compare",
                                        style=style,
                                        value=False,
                                        )

        compare_button.on_click(wrap_compare)
        clear_button = widgets.Button(description="Clear Output",
                                      style=style,
                                      value=False,
                                      )
        clear_button.on_click(lambda a: self.compare(a))
        display(widgets.widgets.HBox(
            (self.molecule_buttons[0], self.molecule_buttons[1])))
        display(widgets.widgets.HBox((compare_button, clear_button)))
