from touscheklib._version import __version__
from touscheklib import optics_tools
import argparse
import json

description = '''touscheklib -- Touschek liftime calculator.

Copyright (C) 2021 Malte Titze
This program comes with ABSOLUTELY NO WARRANTY.
This is free software, and you are welcome to redistribute it
under the GPLv3 license.

'''
epilog=f"v{__version__}"

parser = argparse.ArgumentParser(description=description, epilog=epilog)
parser.add_argument('latticefile', type=str, help="Lattice filename.")
parser.add_argument('beamfile', type=str, help="MAD-X beam input filename in json format.")
parser.add_argument('-p', action='store_true', dest='precise', default=False, help="Precision mode.")
parser.add_argument('--coupling_y', type=float, default=0, dest='coupling_y', help='Coupling factor f'+\
     'to enlarge ey-emittance by ey=f*ex, where ex will be the natural emittance of the lattice.')
parser.add_argument('--symmetry', type=int, default=1, dest='symmetry', help='Machine symmetry.')
parser.add_argument('--n_slices', type=int, default=11, dest='n_slices', help='Number of slices for thick-elements.')
parser.add_argument('--resolution', type=int, default=6, dest='resolution', help='Number of optics function evaluations per [m] in drifts.')
parser_namespace = parser.parse_args()

with open(parser_namespace.beamfile, 'r') as f:
    beam = json.load(f)

# Create optics class with given lattice and beam parameters:
opt = optics_tools.optics(lattice=parser_namespace.latticefile, beam_params=beam, show_init=False, verbose=False)

# Compute Touschek-lifetime:
opt.verbose = True # otherwise no results will be printed to command line
touschek_results = opt.touschek_lifetime(precise=parser_namespace.precise, n_slices=parser_namespace.n_slices,
resolution=parser_namespace.resolution, symmetry=parser_namespace.symmetry, coupling_y=parser_namespace.coupling_y)
