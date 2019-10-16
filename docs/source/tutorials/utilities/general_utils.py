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

# These are some imports/helpers for all of the tutorials
style = json.load(open('./utilities/widget_style.json'))
hyper_dict = json.load(open("./utilities/hyperparameter_presets.json"))
hyper_vals = json.load(open("./utilities/hyperparameter_bounds.json"))


def mask_body_order(hyperparams):
    """ This is for backwards compatibility, as we will soon be moving to using
        "body order" in place of "soap_type"

        Author: R. Cersonsky

        Keywords:
            hyperparams - full or partial dictionary of hyperparameters

        Returns:
            dictionary of hyperparameters suitable for a SOAP object
    """

    soap_types = ["RadialSpectrum", "PowerSpectrum", "BiSpectrum"]
    return {"soap_type": soap_types[hyperparams["body_order"]-1],
            **{h: hyperparams[h] for h in hyperparams if h != 'body_order'}}


def markdown_table_from_dict(d, headers):
    """ Function to generate a markdown table from a dictionary

        Author:  R. Cersonsky

        Keywords:
            d       - dictionary of type {parameter: value}
            headers - optional headers for the table

        Returns:
            markdown string
    """
    return '<table>  <thead><tr>{}</tr></thead><tbody>'.format(''.join(['<th>{}</th>'.format(h) for h in headers]))+''.join(['<tr><td>{}</td>{}</tr>'.format(key.replace('_',' ').title(), ''.join(['<td>{}</td>'.format(round(v,4)) for v in d[key]])) for key in d])+'</tbody></table>'


def compute_kernel(calculator, features1, features2=None, kernel_type='Structure', **kwargs):
    """ Function to calculate similarity kernel

        Author: M. Veit or F. Musil, modified by R. Cersonsky

        Keywords:
            calculator  - soap calculator
            features1   - ASE-Atoms-like object to compute kernel of
            features2   - ASE-Atoms-like object to compute kernel of.
                          Default is features1
            kernel_type - "atom" or "structure"
            kwargs      - hyperparameters of soap calculator

        Returns:
            kernel function
    """
    from rascal.models import Kernel

    kernel = Kernel(representation=calculator, name='Cosine', target_type=kernel_type, zeta=2, **kwargs)
    return kernel(X=features1, Y=features2)

def readme_button():
    """ Function to generate and display a button to show the package README

        Author:  R. Cersonsky

    """
    def show_readme_for_button(b):
        clear_output()
        display(Markdown(str(''.join(open('../README.rst')))))

    button = widgets.Button(description="Show README", style=style)
    output = widgets.Output()
    display(button, output)

    button.on_click(show_readme_for_button)


def _button_template_(options, description, disp_func, default=None):
    """ Function to generate and display a ToggleButton instance for a set of
        options

        Author:  R. Cersonsky

        Keywords:
            options      - list of options
            descriptions - string to describe button
            disp_func    - function object to be called when a button is clicked
            default      - value of the button, needs to be in list of options

        Returns:
            widgets.ToggleButtons Object
    """
    button = widgets.ToggleButtons(
        options=options,
        description=description,
        button_style='',
        style=style,
        value=options[0] if default == None else default
    )
    display(button)
    button.observe(disp_func, 'value')
    return button
