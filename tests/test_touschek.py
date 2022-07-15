from touscheklib import __version__, optics_tools

def test_version():
    assert __version__ == '0.8.0'

def test_cpymad():
    import cpymad.libmadx as l; l.start()

def test_run():
    madx_lattice_filename = 'tests/test_lattice_B2.madx'
    bessy_beam = {'particle': 'ELECTRON',
               'energy': 2.5,
               'ex': 1e-7,
               'ey': 1e-7}
    opt = optics_tools.optics(madx_lattice_filename, beam_params=bessy_beam, 
                              show_init=False, verbose=True)
    opt.get_natural_parameters(n_slices=32, resolution=64)
