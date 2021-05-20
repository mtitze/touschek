from touschek import __version__


def test_version():
    assert __version__ == '0.1.0'

def test_cpymad():
    import cpymad.libmadx as l; l.start()

def test_run():
    from touschek import touschek as ts

    madx_lattice_filename = 'tests/test_lattice1.madx'
    bessy3_beam = {'particle': 'ELECTRON',
               'energy': 2.5,
               'ex': 1e-7,
               'ey': 1e-7}

    m = ts.init_madx(lattice=madx_lattice_filename, **bessy3_beam)

    twiss = m.twiss(betx=1, bety=1, alfx=0.1, alfy=0.1)
