# Script to compute the Touschek lifetime.

import numpy as np
import argparse
from scipy import integrate, constants
from scipy.special import iv
import progressbar
import mpmath as mp

from touschek import dee_to_dpp

'''
References:
[1] A. Piwinski: "THE TOUSCHEK EFFECT IN STRONG FOCUSING STORAGE RINGS", DESY 98-179 (1998).
'''

if __name__ == "__main__":
    from _version import __version__
    description = 'Touschek liftime calculator', 
    epilog=f"v{__version__} by Malte Titze"
    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument('latticefile', type=str, help="lattice filename")
    parser.add_argument('--beamfile', dest='beamfile', default='', type=str, help="MAD-X beam input filename")
    parser_namespace = parser.parse_args()
    lattice_filename = parser_namespace.latticefile
    # need to initiate optics class here

def F_integrand(kappa, kappa_m, b1, b2, kappatol=1e-5, exptol=500):
    '''
    The integrand F in Piwinski's paper [1], Eq. (42)f. This routine is using tolerances in trying to avoid
    NaNs caused e.g. by large internal exponents and may be inaccurate in some circumstances. For counter-check, 
    please use `F_integrand_mp`.

    Note that kappa must be smaller than pi/2. The limits are:
    F_integrand(kappa = np.inf) = 0 and 
    F_integrand(kappa = 0) = -np.inf

    Parameters
    ----------
    kappa: float
        The kappa parameter in [1].
    kappa_m: float
        The lower bound of the kappa parameter in [1].
    b1: float
        B1 value in [1].
    b2: float
        B2 value in [1].
    kappatol: float, optional
        Tolerance below pi/2 in computing tau = tan(kappa).
    exptol: float, optional
        Tolerance to set exp(-b1*tau)*J0(b2*tau) to zero if b1*tau and b2*tau are both larger than exptol.
    '''
    if kappa > np.pi/2 - kappatol:
        # prevent infinity at np.pi/2
        return 0
    tau = np.tan(kappa)**2
    tau_m = np.tan(kappa_m)**2
    term1 = ((2 + tau)**2*(tau/tau_m/(1 + tau) - 1)/tau + tau - \
        np.sqrt(tau*tau_m*(1 + tau)) - (2 + 1/2/tau)*np.log(tau/tau_m/(1 + tau)))*np.sqrt(1 + tau)
    if b1*tau > exptol and b2*tau > exptol:
        # prevent inf if inserting too large values in the Bessel function or the exponent.
        return 0
    term2 = term1*np.exp(-b1*tau)
    if term2 == 0:
        return 0
    else:
        return term2*iv(0, b2*tau)

def F_integrand_mp(kappa, kappa_m, b1, b2):
    '''
    mpmath-version of the routine `F_integrand`, now using arbitrary number precision. This routine should 
    be more accurate but will in general be slower than `F_integrand`.

    For further documentation see `F_integrand` routine.
    '''
    tau = mp.tan(kappa)**2
    tau_m = mp.tan(kappa_m)**2
    return ((2 + tau)**2*(tau/tau_m/(1 + tau) - 1)/tau + tau - \
        mp.sqrt(tau*tau_m*(1 + tau)) - (2 + 1/2/tau)*mp.log(tau/tau_m/(1 + tau)))*mp.sqrt(1 + tau)*mp.exp(-b1*tau)*mp.besseli(0, b2*tau)

def prepare_touschek(optics, delta_pm, verbose=True):
    '''
    Compute necessary parameters for Touschek-lifetime calculation formula (42)f in Ref. [1].

    Parameters
    ----------
    optics: :obj: optics
        Optics class containing the lattice and utilities.
    delta_pm: float
        delta_pm parameter in [1], representing the momentum acceptance of the ring.
    verbose: bool, optional
        Verbose mode.

    Returns
    -------
    dict
        Dictionary containing the values B1, B2 as well as the constant in front of the Touschek-lifetime
        equation, kappa_m and the position.
    '''
    gamma0 = optics.beam.gamma.value
    beta0 = optics.beam.beta.value
    N_p = optics.beam.npart.value
    epsilon_x = optics.beam.ex.value
    epsilon_y = optics.beam.ey.value
    delta_p = dee_to_dpp(optics.beam.sige.value, beta0=beta0)
    sigma_s = optics.beam.sigt.value

    tau_m = beta0**2*delta_pm**2
    kappa_m = np.arctan(np.sqrt(tau_m))

    if verbose:
        print ('\n*** Touschek lifetime input parameters ***')
        print (f'     Np: {N_p}')
        print (f' ex_rms: {epsilon_x}')
        print (f' ey_rms: {epsilon_y}')
        print (f'   dp/p: {delta_p}')
        print (f'sigma_s: {sigma_s}\n')
        print ('*** Optics resolution ***')
        print (f"        n_slices: {optics.function_parameters['n_slices']}")
        print (f"drift resolution: {optics.function_parameters['resolution']}")

    rp = constants.physical_constants['classical electron radius'][0]

    alphax = optics.function['alphax'].values
    alphay = optics.function['alphay'].values
    betax = optics.function['betax'].values
    betay = optics.function['betay'].values
    dx = optics.function['dispx'].values
    dy = optics.function['dispy'].values
    dpx = optics.function['ddispx'].values
    dpy = optics.function['ddispy'].values
    position = optics.function['position'].values

    tilde_dx = alphax*dx + betax*dpx
    tilde_dy = alphay*dy + betay*dpy

    sigma_betax = np.sqrt(betax*epsilon_x + dx**2*delta_p**2) # sigma_x**2 = beta_x*epsilon_x + D_x**2*delta_p**2
    sigma_betay = np.sqrt(betay*epsilon_y + dy**2*delta_p**2)

    sigma_x = np.sqrt(sigma_betax**2 + delta_p**2*dx**2)
    sigma_y = np.sqrt(sigma_betay**2 + delta_p**2*dy**2)

    #tilde_sigma_x2 = sigma_betax**2 + delta_p**2*(dx**2 + tilde_dx**2)
    #tilde_sigma_y2 = sigma_betay**2 + delta_p**2*(dx**2 + tilde_dy**2)

    sigma_h2 = 1/(1/delta_p**2 + (dx**2 + tilde_dx**2)/sigma_betax**2 + (dy**2 + tilde_dy**2)/sigma_betay**2)

    B1 = 1/(2*beta0**2*gamma0**2)*((sigma_betax**2 - sigma_h2*tilde_dx**2)/epsilon_x**2 + \
         (sigma_betay**2 - sigma_h2*tilde_dy**2)/epsilon_y**2)

    B12mB22 = betax**2*betay**2*sigma_h2/(beta0**4*gamma0**4*sigma_betax**4*sigma_betay**4*delta_p**2)*\
                (sigma_x**2*sigma_y**2 - delta_p**4*dx**2*dy**2) # B1**2 - B2**2

    touschek_const = rp**2*constants.c*N_p/(4*np.sqrt(np.pi)*gamma0**4*beta0**2)*\
        np.sqrt(sigma_h2)/(sigma_s*delta_p*epsilon_x*epsilon_y)

    B2 = np.sqrt(B1**2 - B12mB22)
    return {'B1': B1, 'B2': B2, 'const': touschek_const, 'kappa_m': kappa_m,
            'position': position}


def lifetime(precise=False, precision=16, verbose=True, symmetry=1, **kwargs):
    '''
    Compute the Touschek lifetime according to Eq. (42)f, Ref [1].

    Parameters
    ----------
    precise: bool, optional 
        If True, use arbitrary number precision.
    precision: int, optional 
        If precise==True, then this will define the number of digits to be considered.
    symmetry: int, optional
        If the machine admits a known symmetry, one can speed up the calculation process
        by calculating only the first 1/symmetry-th part of the whole machine.
    verbose: bool, optional
        Verbose mode.
    **kwargs
        Optional arguments given to `prepare_touschek`.
    '''

    params = prepare_touschek(verbose=verbose, **kwargs)

    B1 = params['B1']
    B2 = params['B2']
    touschek_const = params['const']
    kappa_m = params['kappa_m']
    position = params['position']

    if verbose and precise:
        print ('\n*** Arbitrary number precision mode ***')
        print (f'Precision: {precision}')


    if symmetry > 1:
        # use only the first part of the symmetric machine
        circumference = position[-1]
        indices = position <= circumference/symmetry
        B1 = B1[indices]
        B2 = B2[indices]
        touschek_const = touschek_const[indices]
        position = position[indices]

    print (flush=True) # clean line before showing progress bar

    touschek_ring = []
    bar = progressbar.ProgressBar(maxval=len(B1), \
        widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
    bar.start()
    if not precise:
        for l in range(len(B1)):
            bar.update(l + 1)
            touschek_ring.append(integrate.quad(lambda x: F_integrand(x, kappa_m=kappa_m, b1=B1[l], b2=B2[l]), kappa_m, np.pi/2)[0])
    else:
        mp.dps = precision
        km = mp.mpf(kappa_m)
        b1m = mp.mpf(1)*B1
        b2m = mp.mpf(1)*B2
        for l in range(len(B1)):
            bar.update(l + 1)
            integrand_mp = lambda x: F_integrand_mp(x, kappa_m=km, b1=b1m[l], b2=b2m[l])
            touschek_ring.append(float(mp.quad(integrand_mp, [km, mp.pi/2])))
    bar.finish()

    ds = np.diff(position)
    touschek_ring = np.array(touschek_ring)
    # n.b. the average will be sufficient also for just a part of a symmetric machine:
    touschek_lifetime_i = sum(touschek_const[1:]*touschek_ring[1:]*ds)/sum(ds)

    if symmetry > 1:
        # combine the result back to the whole machine
        B1 = np.array(list(B1)*symmetry)
        B2 = np.array(list(B2)*symmetry)
        touschek_ring = np.array(list(touschek_ring)*symmetry)
        touschek_const = np.array(list(touschek_const)*symmetry)
        position = np.hstack([position + k for k in np.arange(0, circumference, position[-1])[:-1]])

    lifetime = 1/touschek_lifetime_i

    if verbose:
        print (f'Touschek lifetime [s]: {lifetime:.3f}')
        print (f'                  [h]: {lifetime/60/60:.3f}')

    return {'lifetime': lifetime, 'touschek_ring': touschek_ring, 'symmetry': symmetry,
    'B1': B1, 'B2': B2, 's': position, 'touschek_const': touschek_const, 'kappa_m': kappa_m}
    