from touschek._version import __version__
import numpy as np

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

#import argparse
#if __name__ == "__main__":
#    description = 'Touschek liftime calculator', 
#    epilog=f"v{__version__} by Malte Titze"
#    parser = argparse.ArgumentParser(description=description, epilog=epilog)
#    parser.add_argument('latticefile', type=str, help="lattice filename")
#    parser.add_argument('--beamfile', dest='beamfile', default='', type=str, help="MAD-X beam input filename")
#    parser_namespace = parser.parse_args()
#    lattice_filename = parser_namespace.latticefile
#    # need to initiate optics class here
