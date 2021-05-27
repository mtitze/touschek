# Script to compute the touschek lifetime.

import numpy as np


if __name__ == "__main__":
    from prepare import parse, get_optics_functions, init_madx
    from _version import __version__
    description = 'Touschek liftime calculator', 
    epilog=f"v{__version__} by Malte Titze"
    parser = parse(description=description, epilog=epilog)
    parser_namespace = parser.parse_args()
    lattice_filename = parser_namespace.latticefile




lattice_filename = 'tests/bessy3_5ba-20p_v_long-bend-tgrb.madx'
# run this in interactive jupyter for test:
#lattice_filename = '../tests/bessy3_5ba-20p_v_long-bend-tgrb.madx'


bessy3_beam = {'particle': 'ELECTRON',
               'energy': 2.5,
               'ex': 1e-7,
               'ey': 1e-7}

madx = init_madx(lattice_filename, **bessy3_beam)

optics = get_optics_functions(madx=madx, n_slices=8, resolution=21)

print (madx.table.summ.q1)
#print (m.table.ptc_twiss_summary.q1)
print ('-----')

s = optics['pos']
betax = optics['betax']
ds = np.diff(s)
qx_twiss = sum(1/betax[1:]*ds)/(2*np.pi)
print (qx_twiss)

