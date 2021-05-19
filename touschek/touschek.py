# Script to compute the touschek lifetime.

def init_madx(lattice, **beam_params):
    '''
    Create cpymad instance with given lattice and beam parameters.
    '''
    from cpymad.madx import Madx
    m = Madx()
    m.beam(**beam_params)
    m.input(f'CALL, FILE="{lattice}";')
    return m

if __name__ == "__main__":
    from prepare import parse
    from _version import __version__
    description = 'Touschek liftime calculator', 
    epilog=f"v{__version__} by Malte Titze"
    parse(description=description, epilog=epilog)

# example; move lattice and this to test later on ... 
madx_lattice_filename = 'tests/bessy3_5ba-20p_v_long-bend-tgrb.madx'

bessy3_beam = {'particle': 'ELECTRON',
               'energy': 2.5,
               'ex': 1e-7,
               'ey': 1e-7}

m = init_madx(madx_lattice_filename, **bessy3_beam)

twiss = m.twiss()
