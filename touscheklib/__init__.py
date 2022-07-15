import numpy as np
import toml
import os

tomldict = toml.load(f'{os.path.dirname(__file__)}/../pyproject.toml')
__version__ = tomldict['tool']['poetry']['version']

def dee_to_dpp(dee, beta0):
    '''
    Convert dE/E to dp/p in the absence of an electric potential.
    '''
    return np.sqrt((dee + 1)**2/beta0**2 + 1.0 - 1/beta0**2) - 1

def dpp_to_dee(dpp, beta0):
    '''
    Convert dp/p to dE/E in the absence of an electric potential.
    '''
    return np.sqrt((dpp + 1)**2 + 1/beta0**2 - 1)*beta0 - 1
