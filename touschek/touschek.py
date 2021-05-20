# Script to compute the touschek lifetime.

# %%

def init_madx(lattice, **beam_params):
    '''
    Create cpymad instance with given lattice and beam parameters.
    '''
    from cpymad.madx import Madx
    m = Madx()
    m.beam(**beam_params)
    m.call(lattice)
    return m

# %%

if __name__ == "__main__":
    from prepare import parse
    from _version import __version__
    description = 'Touschek liftime calculator', 
    epilog=f"v{__version__} by Malte Titze"
    parser = parse(description=description, epilog=epilog)
    parser_namespace = parser.parse_args()
    lattice_filename = parser_namespace.latticefile
# %%

#lattice_filename = 'tests/bessy3_5ba-20p_v_long-bend-tgrb.madx'
# run this in interactive jupyter for test:
#lattice_filename = '../tests/bessy3_5ba-20p_v_long-bend-tgrb.madx'


bessy3_beam = {'particle': 'ELECTRON',
               'energy': 2.5,
               'ex': 1e-7,
               'ey': 1e-7}

m = init_madx(lattice_filename, **bessy3_beam)

twiss = m.twiss()# betx=1, bety=1, alfx=0.1, alfy=0.1)

# see http://mad.web.cern.ch/mad/madx.old/ptc_general/ptc_general.html for details.
#ptc_nst = 4
ptc_nst = 50

m.ptc_create_universe()
m.ptc_create_layout(time=True, model=1, method=2, nst=ptc_nst, exact=True)
m.select(flag='PTC_TWISS', clear=True)

# see http://mad.web.cern.ch/mad/madx.old/ptc_twiss/ptc_twiss.html
m.ptc_twiss(icase=6, no=2, deltap=0.0, closed_orbit=True, deltap_dependency=True,
            slice_magnets=True, file='ptc_twiss')
m.ptc_end()