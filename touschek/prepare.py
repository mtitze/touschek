# Script to collect MAD-X beam and MAD-X lattice input.

def parse(**kwargs):
    import argparse
    parser = argparse.ArgumentParser(**kwargs)
    parser.add_argument('latticefile', type=str, help="lattice filename")
    parser.add_argument('--beamfile', dest='beamfile', default='', type=str, help="MAD-X beam input filename")
    return parser

def init_madx(lattice: str, show_init=True, stdout=True, **beam_params):
    '''
    Create cpymad instance with given lattice and beam parameters.
    '''
    from cpymad.madx import Madx
    m = Madx(stdout=stdout)
    if not show_init:
        m.option(echo=False, warn=False)
    m.beam(**beam_params)
    m.call(lattice)
    m.option(echo=True, warn=True)
    return m


import numpy as np

def beta_drift(s, alpha, gamma):
    '''
    Compute the beta-function of a drift if alpha and gamma are known at position s.
    Also return s_star and gamma_star.
    '''
    beta_star = 1/gamma
    s_star = s + beta_star*alpha
    return lambda x: beta_star + (x - s_star)**2/beta_star, s_star, beta_star

def dispersion_drift(s, disp, ddisp):
    '''
    Return the dispersion function and d(dispersion)/ds-function of a drift.
    '''
    return lambda x: ddisp*(x - s) + disp, lambda x: ddisp*np.ones(len(x))

def insert_optics_values(twiss_table, resolution):
    '''
    Compute missing beta-values of the drift section of a MAD-X twiss table.
    '''
    
    drift_elements=['drift', 'marker', 'instrument', 'monitor', 'monitorv', 'monitorh']
    keywords = twiss_table.keyword
    lengths = twiss_table.L
    positions = twiss_table.s

    drift_indices = np.array([k for k in range(len(keywords)) if keywords[k] in drift_elements and lengths[k] > 0])
    
    additional_positions = []
    additional_betax, additional_betay = [], []
    additional_alphax, additional_alphay = [], []
    additional_gammax, additional_gammay = [], []
    additional_dx, additional_dpx = [], []
    additional_dy, additional_dpy = [], []

    for k in drift_indices:
        n_points = int(resolution*lengths[k])
        base = np.linspace(positions[k] - lengths[k]/2, 
                           positions[k] + lengths[k]/2, 
                           n_points)[1:-1]
        additional_positions.append(base)

        betax_func, sx_star, betax_star = beta_drift(positions[k], alpha=twiss_table.alfa11[k], gamma=twiss_table.gama11[k])
        additional_betax.append(betax_func(base))

        betay_func, sy_star, betay_star = beta_drift(positions[k], alpha=twiss_table.alfa22[k], gamma=twiss_table.gama22[k])
        additional_betay.append(betay_func(base))

        # In drift it holds: gamma = const. =: gamma_0 and therefore alpha = -gamma_0*(s - s_0) + alpha_0.
        additional_gammax.append([twiss_table.gama11[k]]*len(base))
        additional_gammay.append([twiss_table.gama22[k]]*len(base))
        additional_alphax.append(twiss_table.gama11[k]*(positions[k] - base) + twiss_table.alfa11[k])
        additional_alphay.append(twiss_table.gama22[k]*(positions[k] - base) + twiss_table.alfa22[k])

        dx_func, dpx_func = dispersion_drift(positions[k], disp=twiss_table.dx[k], ddisp=twiss_table.dpx[k])
        dy_func, dpy_func = dispersion_drift(positions[k], disp=twiss_table.dy[k], ddisp=twiss_table.dpy[k])

        additional_dx.append(dx_func(base))
        additional_dpx.append(dpx_func(base))
        additional_dy.append(dy_func(base))
        additional_dpy.append(dpy_func(base))

    new_positions = np.concatenate((positions, np.hstack(additional_positions)))
    sindices = np.argsort(new_positions)

    out = {}
    out['position'] = new_positions[sindices]
    out['betax'] = np.concatenate((twiss_table.beta11, np.hstack(additional_betax)))[sindices]
    out['betay'] = np.concatenate((twiss_table.beta22, np.hstack(additional_betay)))[sindices]
    out['alphax'] = np.concatenate((twiss_table.alfa11, np.hstack(additional_alphax)))[sindices]
    out['alphay'] = np.concatenate((twiss_table.alfa22, np.hstack(additional_alphay)))[sindices]
    out['gammax'] = np.concatenate((twiss_table.gama11, np.hstack(additional_gammax)))[sindices]
    out['gammay'] = np.concatenate((twiss_table.gama22, np.hstack(additional_gammay)))[sindices]

    out['dispx'] = np.concatenate((twiss_table.dx, np.hstack(additional_dx)))[sindices]
    out['ddispx'] = np.concatenate((twiss_table.dpx, np.hstack(additional_dpx)))[sindices]
    out['dispy'] = np.concatenate((twiss_table.dy, np.hstack(additional_dy)))[sindices]
    out['ddispy'] = np.concatenate((twiss_table.dpy, np.hstack(additional_dpy)))[sindices]
    return out

def get_optics_functions(madx, n_slices=41, resolution=201):
    '''
    Compute the optics functions for a given lattice, using the MAD-X twiss functionality.
    
    INPUT
    =====
    n_slices: number of slices to split elements when computing the MAD-X thin lattice
    resolution: number of points/m for drift sections.
    '''

    madx.select(flag='makethin', clear=True)
    madx.select(flag='makethin', thick=True, slice=n_slices)
    madx.makethin(style='simple', sequence='ring')
    madx.use(sequence='ring')
    thin_twiss = madx.twiss(centre=True, ripken=True, chrom=True)

    return insert_optics_values(thin_twiss, resolution=resolution)