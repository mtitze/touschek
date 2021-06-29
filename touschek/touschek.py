# Script to compute the Touschek lifetime.

import numpy as np
import argparse
from scipy import integrate, constants
from scipy.special import iv
import progressbar
import mpmath as mp

from touschek import dee_to_dpp

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

def F_integrand(kappa, kappa_m, b1, b2):
    '''
    Note that kappa must be smaller than pi/2. The limits are:
    F_integrand(kappa = np.inf) = 0 and 
    F_integrand(kappa = 0) = -np.inf
    '''
    tol = 1e-5
    if kappa > np.pi/2 - tol:
        # prevent infinity at np.pi/2
        return 0
    tau = np.tan(kappa)**2
    tau_m = np.tan(kappa_m)**2
    term1 = ((2 + tau)**2*(tau/tau_m/(1 + tau) - 1)/tau + tau - \
        np.sqrt(tau*tau_m*(1 + tau)) - (2 + 1/2/tau)*np.log(tau/tau_m/(1 + tau)))*np.sqrt(1 + tau)
    if b1*tau > 500 and b2*tau > 500:
        # prevent inf if inserting too large values in the Bessel function or the exponent.
        return 0
    term2 = term1*np.exp(-b1*tau)
    if term2 == 0:
        return 0
    else:
        return term2*iv(0, b2*tau)

def F_integrand_mp(kappa, kappa_m, b1, b2):
    '''
    mpmath-version of the Touschek integrand for arbitrary number precision.

    Note that kappa must be smaller than pi/2. The limits are:
    F_integrand(kappa = np.inf) = 0 and 
    F_integrand(kappa = 0) = -np.inf
    '''
    tau = mp.tan(kappa)**2
    tau_m = mp.tan(kappa_m)**2
    return ((2 + tau)**2*(tau/tau_m/(1 + tau) - 1)/tau + tau - \
        mp.sqrt(tau*tau_m*(1 + tau)) - (2 + 1/2/tau)*mp.log(tau/tau_m/(1 + tau)))*mp.sqrt(1 + tau)*mp.exp(-b1*tau)*mp.besseli(0, b2*tau)

def prepare_touschek(optics, delta_pm, verbose=True):
    '''
    Compute the necessary parameters for Touschek-lifetime calculation.
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


def lifetime(precise=False, precision=16, verbose=True, **kwargs):
    '''
    Compute the Touschek lifetime.

    precise: If True, use arbitrary number precision to deal with the integrand.
    precision: If precise==True, then this will denote the number of digits to be considered.
    '''

    params = prepare_touschek(verbose=verbose, **kwargs)

    B1 = params['B1']
    B2 = params['B2']
    touschek_const = params['const']
    kappa_m = params['kappa_m']

    if verbose and precise:
        print ('\n*** Arbitrary number precision mode ***')
        print (f'Precision: {precision}')

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

    position = params['position']
    ds = np.diff(position)
    touschek_ring = np.array(touschek_ring)
    touschek_lifetime_i = sum(touschek_const[1:]*touschek_ring[1:]*ds)/sum(ds)
    lifetime = 1/touschek_lifetime_i

    if verbose:
        print (f'Touschek lifetime [s]: {lifetime:.3f}')
        print (f'                  [h]: {lifetime/60/60:.3f}')

    return {'lifetime': lifetime, 'touschek_ring': touschek_ring, 
    'B1': B1, 'B2': B2, 's': position, 'touschek_const': touschek_const, 'kappa_m': kappa_m}
    