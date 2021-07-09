#from typing import Sequence
from os import sync
from cpymad.madx import Madx # BaseTypeMap
import numpy as np
import pandas as pd
from scipy import constants
import itertools

import warnings
import touschek

from touschek.plotting import plot_survey, plot_touschek_losses
from touschek.touschek import lifetime
from touschek import dee_to_dpp

def init_madx(lattice: str, show_init=True, verbose=True, **kwargs):
    '''
    Create cpymad instance with given lattice and beam parameters.
    '''
    if verbose:
        print ('loading lattice:\n{}'.format(lattice))
    m = Madx(**kwargs)
    if not show_init:
        m.option(echo=False, warn=False)
    m.call(lattice)
    m.option(echo=True, warn=True)
    return m

def get_beam_parameters(madx):
    '''
    Get current madx beam parameters.
    '''
    #Obtain beam parameters bye defining global variables in order 
    #to be able to load them from the madx instance.
    madx.option(info=False)
    madx.input('beam_energy := beam->energy;')
    madx.input('beam_gamma := beam->gamma;')
    madx.input('beam_beta := beam->beta;')
    madx.input('beam_particle_mass := beam->mass;')
    madx.input('beam_charge := beam->charge;')
    madx.input('beam_ex := beam->ex;')
    madx.input('beam_ey := beam->ey;')
    madx.input('beam_p0c := beam->pc;')
    madx.input('beam_exn := beam->exn;')
    madx.input('beam_eyn := beam->eyn;')
    madx.input('beam_npart := beam->npart;')
    madx.input('beam_sige := beam->sige;')
    madx.input('beam_sigt := beam->sigt;')
    #madx.input('beam_particle := beam->particle;')
    madx.option(info=True)
    
    return {'energy': {'value': madx.globals.beam_energy, 'description': 'particle energy in GeV'},
            'gamma': {'value': madx.globals.beam_gamma, 'description': 'relativistic gamma factor'},
            'beta': {'value': madx.globals.beam_beta, 'description': 'relativistic beta factor'},
            'mass': {'value': madx.globals.beam_particle_mass, 'description': 'particle mass in GeV'},
            'charge': {'value': madx.globals.beam_charge, 'description': 'particle charge in [e]'},
            'pc': {'value': madx.globals.beam_p0c, 'description': 'p0*c: particle momentum * speed of light'},
            'ex': {'value': madx.globals.beam_ex, 'description': 'rms x-emittance [m]'},
            'ey': {'value': madx.globals.beam_ey, 'description': 'rms y-emittance [m]'},
            'exn': {'value': madx.globals.beam_exn, 'description': 'energy-normalized rms x-emittance [m]'},
            'eyn': {'value': madx.globals.beam_eyn, 'description': 'energy-normalized rms y-emittance [m]'},
            'npart': {'value': madx.globals.beam_npart, 'description': 'number of particles'},
            'sige': {'value': madx.globals.beam_sige, 'description': 'relative energy spread sigma_E/E'},
            'sigt': {'value': madx.globals.beam_sigt, 'description': 'rms bunch length [m]'}
            }
            #'particle': {'value': madx.beam.particle, 'description': 'particle type'}

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

def test_tune(optics):
    print ('TUNE TEST')
    print ('---------')
    q1 = optics.thin_twiss.summary.q1
    q2 = optics.thin_twiss.summary.q2
    print (f'qx-goal: {q1}')
    print (f'qy-goal: {q2}')
    print ('---------')
    s = optics.function['position'].values
    betax = optics.function['betax'].values
    betay = optics.function['betay'].values
    ds = np.diff(s)
    qx = sum(1/betax[1:]*ds)/(2*np.pi)
    qy = sum(1/betay[1:]*ds)/(2*np.pi)
    print (f'qx-result: {qx}')
    print (f'qy-result: {qy}')
    print ('Differences:')
    print (f'|qx-goal - qx-result| = {abs(q1 - qx)}')
    print (f'|qy-goal - qy-result| = {abs(q2 - qy)}')


from dataclasses import dataclass
@dataclass
class beam_parameter:
    '''
    Generic class to hold a beam parameter value
    '''
    key: str
    value: float 
    description: str = ''

class beam:
    '''
    Namespace to hold the beam parameters
    '''
    def __init__(self, description=''):
        self.description = description

class optics:
    def __init__(self, lattice: str, beam_params: dict, sequence_name='ring', verbose=True, **kwargs):
        self.lattice = lattice
        self.sequence_name = sequence_name
        self.verbose = verbose
        self.madx = init_madx(lattice=lattice, verbose=verbose, **kwargs)
        self.madx.beam(**beam_params)
        self.get_beam_parameters(madx=self.madx)

    def get_beam_parameters(self, madx=None, update_only_if_different=True):
        '''
        Get beam parameters from given MAD-X instance.

        update_only_if_different: If True, then update existing beam parameters only if they differ
        from the current MAD-X instance.
        '''
        if madx == None:
            madx = self.madx

        madx_beam_parameters = get_beam_parameters(madx=madx)

        if not hasattr(self, 'beam'):
            if self.verbose:
                print (f"Creating beam parameter namespace.")
            setattr(self, 'beam', beam())
        else:
            if self.verbose:
                print (f"Existing beam parameter namespace found.")

        for key in madx_beam_parameters:
            value = madx_beam_parameters[key]['value']
            if hasattr(self.beam, key):
                old_value = getattr(self.beam, key).value
                if update_only_if_different and old_value == value:
                    continue
                if self.verbose:
                    print ('Updating parameter {} ... OLD: {} NEW: {}'.format(key, old_value, value))
            bp = beam_parameter(key=key, value=value, description=madx_beam_parameters[key]['description'])
            setattr(self.beam, key, bp)

    def set_beam_parameters(self, madx=None):
        '''
        Set beam parameters to given MAD-X instance.
        '''
        if madx == None:
            madx = self.madx
        if not hasattr(self, 'beam'):
            raise ValueError(f"ERROR: beam parameter namespace not found.")

        madx_beam_params = {}
        for key in ['energy', 'gamma', 'beta', 'mass', 'charge', 'pc', 
                    'ex', 'ey', 'exn', 'eyn', 'npart', 'sige','sigt']:
            value = getattr(self.beam, key).value
            madx_beam_params[key] = value

        madx.beam(**madx_beam_params)

    def twiss(self, madx=None, centre=True, ripken=True, chrom=True, rmatrix=True, **kwargs):
        '''
        Invoke MAD-X twiss with standard parameters of our interest.
        '''
        if madx == None:
            madx = self.madx
        madx.use(sequence=self.sequence_name)
        madx.select(flag='twiss', clear=True)
        return madx.twiss(centre=centre, ripken=ripken, chrom=chrom, rmatrix=rmatrix, **kwargs)

    def ptc_twiss(self, madx=None, time=True, model=1, method=2, nst=6, exact=True,
                  icase=6, no=2, deltap=0.0, closed_orbit=True,
                  deltap_dependency=True, ring_parameters=True, maptable=True):
        '''
        Invoke MAD-X PTC_TWISS with standard parameters of our interest.

        Note that chromaticities are computed only for icase = 5.
        '''

        if madx == None:
            madx = self.madx
        madx.use(sequence=self.sequence_name)

        madx.ptc_create_universe()
        madx.ptc_create_layout(time=time, model=model, method=method, nst=nst, exact=exact)
        madx.select(flag='PTC_TWISS', clear=True)
        # see http://mad.web.cern.ch/mad/madx.old/ptc_twiss/ptc_twiss.html
        madx.ptc_twiss(icase=icase, no=no, deltap=deltap, 
        closed_orbit=closed_orbit, deltap_dependency=deltap_dependency,
                        ring_parameters=ring_parameters, maptable=maptable,
                        table='ptc_twiss_table')#, file='ptc_twiss')
        madx.ptc_end()
        return madx.table.ptc_twiss_table

    def makethin(self, n_slices, style='simple', stdout=False, **kwargs):
        '''
        MAD-X makethin commands. To prevent interfering with the original lattice, a new MAD-X instance
        is created.
        (internal MAD-X commands like 'extract' and 'use' are unreliable. E.g. 'use' of original lattice 
        after 'makethin' on an extracted sequence does not work) .
        '''
        if self.verbose:
            print ('Creating separate MAD-X instance for makethin ...')
        self.madx_thin = init_madx(lattice=self.lattice, show_init=False, stdout=stdout)
        self.set_beam_parameters(madx=self.madx_thin)
        self.madx_thin.use(sequence=self.sequence_name)
        self.madx_thin.select(flag='makethin', clear=True)
        self.madx_thin.select(flag='makethin', thick=True, slice=n_slices)
        self.madx_thin.makethin(style=style, sequence=self.sequence_name, **kwargs)
        self.madx_thin.use(sequence=self.sequence_name)
        if self.verbose:
            print ('done.')

    def get_one_turn_maps(self, twiss_table=None):
        '''
        Return list of one-turn-maps at each twiss position.
        '''
        if twiss_table == None and hasattr(self.madx.table, 'twiss'):
            twiss_table = self.madx.table.twiss
        else:
            raise ValueError('Twiss table not found.')

        twiss_positions = len(twiss_table.s)
        one_turn_maps = np.zeros([6, 6, twiss_positions])
        for i, j in itertools.product(np.arange(6), repeat=2):
            one_turn_maps[i, j, :] = getattr(twiss_table, f're{i + 1}{j + 1}')
        return one_turn_maps

    def compute_optics_functions(self, n_slices=6, resolution=11, style='simple', **kwargs):
        '''
        Compute the optics functions for a given lattice, using the MAD-X twiss functionality.
        
        INPUT
        =====
        n_slices: number of slices to split elements when computing the MAD-X thin lattice
        resolution: number of points/m for drift sections.

        Recommendation:
        For quick tests, n_slices=6, resolution=11 should be fine. If more precision is
        required, for example to determine the tune, then n_slices=41, resolution=201 turned
        out to be sufficient.
        '''
        if self.verbose:
            print ('Computing optics functions with the following parameters:')
            print (f'  n_slices: {n_slices}')
            print (f'resolution: {resolution}')

        self.makethin(n_slices=n_slices, style=style, **kwargs)

        if self.verbose:
            print ('Performing MAD-X twiss on thin lattice ...')
        self.thin_twiss = self.twiss(madx=self.madx_thin)

        keywords = self.thin_twiss.keyword
        lengths = self.thin_twiss.L
        positions = self.thin_twiss.s

        drift_elements=['drift', 'marker', 'instrument', 'monitor', 'monitorv', 'monitorh']
        drift_indices = np.array([k for k in range(len(keywords)) if keywords[k] in drift_elements and lengths[k] > 0])
        
        # also compute 1/rho
        inv_radii_x = np.zeros(len(positions))
        inv_radii_y = np.zeros(len(positions))

        non_zero_lengths = np.where(np.logical_and(lengths != 0, [keyword in ['sbend', 'rbend', 'dipedge'] for
        keyword in keywords]))[0]
        inv_radii_x[non_zero_lengths] = self.thin_twiss.angle[non_zero_lengths]/lengths[non_zero_lengths]
        inv_radii_y[non_zero_lengths] = self.thin_twiss.k0sl[non_zero_lengths]/lengths[non_zero_lengths] # n.b. the vertical angle in MAD-X may be hidden in rotations before and after bends

        # For the fourth synchrotron integral it is required (see Wolski, p. 225) to also use any normal-quad components
        # found inside bends (if they exist, as in combined-function magnets)
        # n.b. k1l in MAD-X twiss correspnds to 1/(B rho)*\partial B_y/\partial x, see Manual p. 75.
        k1bends_x = np.zeros(len(positions))
        k1bends_x[non_zero_lengths] = self.thin_twiss.k1l[non_zero_lengths]
        k1bends_y = np.zeros(len(positions))
        k1bends_y[non_zero_lengths] = self.thin_twiss.k1sl[non_zero_lengths]

        additional_positions = []
        additional_betax, additional_betay = [], []
        additional_alphax, additional_alphay = [], []
        additional_gammax, additional_gammay = [], []
        additional_dx, additional_dpx = [], []
        additional_dy, additional_dpy = [], []
        additional_inv_radii_x = []
        additional_inv_radii_y = []
        additional_k1bends_x = []
        additional_k1bends_y = []

        for k in drift_indices:
            n_points = int(resolution*lengths[k])
            base = np.linspace(positions[k] - lengths[k]/2, 
                            positions[k] + lengths[k]/2, 
                            n_points)[1:-1]
            additional_positions.append(base)

            betax_func, sx_star, betax_star = beta_drift(positions[k], alpha=self.thin_twiss.alfa11[k], gamma=self.thin_twiss.gama11[k])
            additional_betax.append(betax_func(base))

            betay_func, sy_star, betay_star = beta_drift(positions[k], alpha=self.thin_twiss.alfa22[k], gamma=self.thin_twiss.gama22[k])
            additional_betay.append(betay_func(base))

            # In drift it holds: gamma = const. =: gamma_0 and therefore alpha = -gamma_0*(s - s_0) + alpha_0.
            additional_gammax.append([self.thin_twiss.gama11[k]]*len(base))
            additional_gammay.append([self.thin_twiss.gama22[k]]*len(base))
            additional_alphax.append(self.thin_twiss.gama11[k]*(positions[k] - base) + self.thin_twiss.alfa11[k])
            additional_alphay.append(self.thin_twiss.gama22[k]*(positions[k] - base) + self.thin_twiss.alfa22[k])

            dx_func, dpx_func = dispersion_drift(positions[k], disp=self.thin_twiss.dx[k], ddisp=self.thin_twiss.dpx[k])
            dy_func, dpy_func = dispersion_drift(positions[k], disp=self.thin_twiss.dy[k], ddisp=self.thin_twiss.dpy[k])

            additional_dx.append(dx_func(base))
            additional_dpx.append(dpx_func(base))
            additional_dy.append(dy_func(base))
            additional_dpy.append(dpy_func(base))
            
            additional_inv_radii_x.append(np.zeros(len(base)))
            additional_k1bends_x.append(np.zeros(len(base)))
            additional_inv_radii_y.append(np.zeros(len(base)))
            additional_k1bends_y.append(np.zeros(len(base)))

        new_positions = np.concatenate((positions, np.hstack(additional_positions)))
        sindices = np.argsort(new_positions)

        function = pd.DataFrame(columns=['position', 'betax', 'betay', 'alphax', 'alphay', 'gammax', 'gammay',
                                          'dispx', 'ddispx', 'dispy', 'ddispy', 'inv_radius'])

        function['position'] = new_positions[sindices]
        function['betax'] = np.concatenate((self.thin_twiss.beta11, np.hstack(additional_betax)))[sindices]
        function['betay'] = np.concatenate((self.thin_twiss.beta22, np.hstack(additional_betay)))[sindices]
        function['alphax'] = np.concatenate((self.thin_twiss.alfa11, np.hstack(additional_alphax)))[sindices]
        function['alphay'] = np.concatenate((self.thin_twiss.alfa22, np.hstack(additional_alphay)))[sindices]
        function['gammax'] = np.concatenate((self.thin_twiss.gama11, np.hstack(additional_gammax)))[sindices]
        function['gammay'] = np.concatenate((self.thin_twiss.gama22, np.hstack(additional_gammay)))[sindices]
        function['dispx'] = np.concatenate((self.thin_twiss.dx, np.hstack(additional_dx)))[sindices]
        function['ddispx'] = np.concatenate((self.thin_twiss.dpx, np.hstack(additional_dpx)))[sindices]
        function['dispy'] = np.concatenate((self.thin_twiss.dy, np.hstack(additional_dy)))[sindices]
        function['ddispy'] = np.concatenate((self.thin_twiss.dpy, np.hstack(additional_dpy)))[sindices]
        function['inv_radius_x'] = np.concatenate((inv_radii_x, np.hstack(additional_inv_radii_x)))[sindices]
        function['k1bends_x'] = np.concatenate((k1bends_x, np.hstack(additional_k1bends_x)))[sindices]
        function['inv_radius_y'] = np.concatenate((inv_radii_y, np.hstack(additional_inv_radii_y)))[sindices]
        function['k1bends_y'] = np.concatenate((k1bends_y, np.hstack(additional_k1bends_y)))[sindices]

        self.function = function
        self.function_parameters = {'n_slices': n_slices, 'resolution': resolution, 'style': style}

        if self.verbose:
            if any(function['k1bends_x'] != 0) or any(function['k1bends_y'] != 0):
                warnings.warn('Combined-function elements found in the lattice.')

    def get_natural_parameters(self, run_emit=False, **kwargs):
        '''
        Compute natural parameters of the given lattice.

        Run emit: If set to True, update internal MAD-X beam parameters using MAD-X EMIT command.
        '''        
        # run emit to check RF setup
        # get optics functions (required to compute e.g. lower emittance limit and counter-check the
        # synchrotron integrals from MAD-X twiss)
        self.compute_optics_functions(**kwargs) # self.madx_thin created

        if run_emit:
            self.madx.use(sequence=self.sequence_name)
            self.madx.emit() # MAD-X emit has internal modules to obtain the natural rms values which are computed here.
            self.get_beam_parameters()

        s = self.function.position.values
        ds = np.diff(s)
        gamma0 = self.beam.gamma.value
        beta0 = self.beam.beta.value
        rho_inv_x = self.function.inv_radius_x.values
        rho_inv_y = self.function.inv_radius_y.values
        energy = self.beam.energy.value*1e9*constants.elementary_charge # in SI units

        circumference = s[-1]

        # Compute the synchrotron integrals
        Dx = self.function.dispx.values
        Dy = self.function.dispy.values
        Dxp = self.function.ddispx.values
        alphax = self.function.alphax.values
        betax = self.function.betax.values
        betay = self.function.betay.values
        gammax = self.function.gammax.values

        si1 = sum(Dx[1:]*rho_inv_x[1:]*ds)
        si2 = sum(rho_inv_x[1:]**2*ds)
        si3 = sum(abs(rho_inv_x[1:])**3*ds)

        k1bends_x = self.function.k1bends_x.values  # the quadrupole gradient in the dipole fields
        k1bends_y = self.function.k1bends_y.values  # the quadrupole gradient in the dipole fields
        # In MAD-X the quadrupole gradient in the bends needs to be divided by
        # the length ds. Therefore ds canceled for these terms.
        si4 = sum(Dx[1:]*rho_inv_x[1:]*(rho_inv_x[1:]**2*ds + 2*k1bends_x[1:]))
        si4_y = sum(Dy[1:]*rho_inv_y[1:]*(rho_inv_y[1:]**2*ds + 2*k1bends_y[1:]))

        Hx = gammax*Dx**2 + 2*alphax*Dx*Dxp + betax*Dxp**2
        si5 = sum(Hx[1:]*abs(rho_inv_x[1:])**3*ds)

        # DERIVED QUANTITIES
        ####################
        particle_mass = self.beam.mass.value*1e9*constants.e/(constants.speed_of_light**2)
        particle_charge = self.beam.charge.value*abs(constants.elementary_charge)

        Cgamma = particle_charge**2/(3*constants.epsilon_0*particle_mass**4*constants.speed_of_light**8)
        # loss of energy per revolution
        u0 = Cgamma*beta0**3*energy**4*si2/(2*np.pi)
        u0_eV = u0/constants.e

        momentum_compaction = si1/circumference
        gamma_tr = 1/np.sqrt(momentum_compaction) # transition energy
        slip_factor = momentum_compaction - 1/gamma0**2

        # NATURAL VALUES
        ################
        # natural x-emittance, see Wolski p.233 Eq. (7.83)
        Cq = 55/32/np.sqrt(3)*constants.hbar/particle_mass/constants.speed_of_light
        jx = 1 - si4/si2
        natural_emittance_x = Cq*gamma0**2*si5/si2/jx

        # lower limit of y-emittance, see Wolski p234 Eq. (7.85) (from Raubenheimer 1991)
        natural_emittance_y_limit = 13/55*Cq/si2*sum(betay[1:]*np.abs(rho_inv_x[1:])**3*ds)

        # nautral energy spread, see Wolski p.236 Eq. (7.94)
        jz = 2 + si4/si2
        natural_dee = gamma0*np.sqrt(Cq*si3/si2/jz) # ! = deltaE/E_0 see Wiedemann p. 302,
        # and Wolski: E/(p0*c) - 1/beta0 = (E - E0)/(p0*c) = \Delta E/E0*beta0 with E0 = p0*c/beta0
        # therefore:
        natural_dpp = dee_to_dpp(natural_dee, beta0=beta0)

        # in order to get the energy acceptance and the synchrotron tune, we need to load the RF cavity parameters.
        # !!! TODO here we assume a *single* RF cavity
        rf_parameters = self.get_rf_cavity_system()
        n_cavities = len(rf_parameters['voltage'])
        if n_cavities == 0:
            raise RuntimeError('No RF cavity found!')
        rf_lag = np.mean(rf_parameters['phase_rad']) # TODO adjust for sophisticated RF system.
        voltage = sum(rf_parameters['voltage']) # TODO adjust for sophisticated RF system.
        if np.sign(voltage) < 0:
            # if voltage is given with different sign, include this as phase advance and correct
            # voltage accordingly.
            rf_lag += np.pi
            voltage = np.abs(voltage)
        harmonic = np.mean(rf_parameters['harmonic']) # TODO adjust for sophisticated RF system.
        sin_phi = u0_eV/voltage # TODO adjust for sophisticated RF system.

        phi = np.arcsin(sin_phi)
        madx_phi = phi/(2*np.pi)
        energy_acceptance = self.energy_acceptance(phase=rf_lag, harmonic=harmonic, 
                                                   voltage=voltage, slip_factor=slip_factor)
        momentum_acceptance = dee_to_dpp(energy_acceptance, beta0=beta0)

        if self.verbose:
            print ()
            print (f'rf total voltage [MV] = {voltage/(1e6)}')
            print (f'rf cavity lag [rad] = {rf_lag}')
            print (f'phi[rad] = {phi}')
            print (f'phi[degree] = {phi*180/np.pi}')
            print (f'MAD-X LAG correction: phi[rad]/(2*PI) = {madx_phi}')
            print (f'energy acceptance: {energy_acceptance} (assuming a single effective RF cavity)')

        # get the synchrotron frequency (Wiedemann, p. 202)
        t_rev = circumference/(beta0*constants.speed_of_light) # revolution time
        f_rev = 1/t_rev # revolution frequency
        omega_rev = 2*np.pi*f_rev # revolution angular frequency

        synchrotron_tune = np.sqrt(-harmonic*slip_factor*constants.e*voltage*np.cos(rf_lag)/\
        (2*np.pi*beta0**2*energy)) # synchrotron angular frequency; "-" sign because Wiedemann defines the slip
        # factor with a negative sign in comparison to our definition. TODO adjust for sophisticated RF system.
        # n.b. synchrotron_tune = omega_s/omega_rev
        omega_s = synchrotron_tune*omega_rev

        # DAMPING TIMES
        ###############
        si2_y = si2
        jy = 1 - si4_y/si2_y

        # Ref. [Wolski], Eq. (7.64), p. 229
        t_damping_f = 2*energy/u0*t_rev
        t_damping_x = t_damping_f/jx
        t_damping_y = t_damping_f/jy
        t_damping_z = t_damping_f/jz

        if self.verbose:
            print (f'               revolution time [s]: {t_rev}')
            print (f'        revolution frequency [1/s]: {f_rev}')
            print (f'revolution angular frequency [1/s]: {omega_rev}')
            print (f'omega_s = {omega_s} (synchrotron angular frequency)')
            print (f'  |nu_s| = {synchrotron_tune} (synchrotron tune)')
            print (f' # turns = {1/synchrotron_tune} (synchrotron oscillation period)')

            print ('Synchrotron damping times:')
            print (f'damping time x [s]: {t_damping_x}')
            print (f'damping time y [s]: {t_damping_y}')
            print (f'damping time z [s]: {t_damping_z}')
            print ('Robinson damping theorem: 4 = jx + jy + jz: {}'.format(jx + jy + jz))

        # natural bunch length, see Wolski, p.237 Eq. (7.96)
        natural_z = circumference/(2*np.pi)*slip_factor/synchrotron_tune*natural_dpp

        # momentum_compaction = ttab.alfa # = alpha_c (momentum compaction factor)
        #natural_z_v2 = constants.speed_of_light/omega_s*natural_dpp*momentum_compaction # wolski, lecture notes
        #natural_z_v3 = constants.speed_of_light/omega_s*natural_dpp*slip_factor # wiedemann p.303, 
        # but see also Eq. (10.26) p. 376 which has a beta-factor included

        pulse_length_ps = natural_z/constants.speed_of_light*1e12 # the expected pulse length in ps
        #pulse_length_ps_v2 = natural_z_v2/constants.speed_of_light*1e12 # the expected pulse length in ps
        #pulse_length_ps_v3 = natural_z_v3/constants.speed_of_light*1e12 # the expected pulse length in ps

        if self.verbose:
            ttab = self.madx_thin.table.summ

            print ('\nSynchrotron integrals')
            print ('---------------------')
            print (f'I1: {si1}, MAD-X: {ttab.synch_1}')
            print (f'I2: {si2}, MAD-X: {ttab.synch_2}')
            print (f'I3: {si3}, MAD-X: {ttab.synch_3}')
            print (f'I4: {si4}, MAD-X: {ttab.synch_4}')
            print (f'I5: {si5}, MAD-X: {ttab.synch_5}')

            print ('\nDerived quantities')
            print ('------------------')
            print (f'C-gamma [m (GeV)**(-3)]: {Cgamma*(1e9*constants.elementary_charge)**3}')
            print ('Energy loss U0 per turn [keV]:')
            print (f'{u0/constants.elementary_charge*1e-3}, MAD-X: ')
            print (f'Momentum compaction factor: {momentum_compaction}')
            print (f'               Slip factor: {slip_factor}')
            print (f'Transition energy gamma_tr: {gamma_tr}')

            print (f'natural x-emittance: {natural_emittance_x}')
            print (f'natural y-emittance lower bound: {natural_emittance_y_limit}')
            print (f'natural_dp/p: {natural_dpp}')
            print ()
            print (f'natural rms bunch length: {natural_z}')
            #print (f'natural rms bunch length v2: {natural_z_v2}')
            #print (f'natural rms bunch length v3: {natural_z_v3}')

            print (f'expected pulse length [ps]: {pulse_length_ps}, FWHM: {pulse_length_ps*2.355}')
            #print (f'expected pulse length [ps]: {pulse_length_ps_v2}, FWHM: {pulse_length_ps_v2*2.355}')
            #print (f'expected pulse length [ps]: {pulse_length_ps_v3}, FWHM: {pulse_length_ps_v3*2.355}')

        # collect and return results
        return {'circumference': circumference, 'si1': si1, 'si2': si2, 'si3': si3, 'si4': si4, 'si5': si5,
         'Cgamma': Cgamma, 'u0': u0, 'u0_ev': u0_eV, 'momentum_compaction': momentum_compaction,
         'gamma_tr': gamma_tr, 'slip_factor': slip_factor, 'nat_ex': natural_emittance_x,
         'nat_ey_lb': natural_emittance_y_limit, 'nat_dee': natural_dee, 'nat_dpp': natural_dpp,
         'voltage': voltage, 'harmonic': harmonic, 'rf_lag': rf_lag, 'rf_lag_u0': phi, 'rf_lag_u0_madx': madx_phi, 'n_cavities': n_cavities,
         'energy_acceptance': energy_acceptance, 'momentum_acceptance': momentum_acceptance, 't_rev': t_rev, 'f_rev': f_rev, 'omega_rev': omega_rev,
         'omega_s': omega_s, 'synchrotron_tune': synchrotron_tune, 'synchrotron_turns': 1/synchrotron_tune,
         'nat_z': natural_z, 'pulse_length_rms_ps': pulse_length_ps, 'pulse_length_fwhm_ps': pulse_length_ps*2.355,
         't_damping_x': t_damping_x, 't_damping_y': t_damping_y, 't_damping_z': t_damping_z,
         'n_slices': self.function_parameters['n_slices'], 'resolution': self.function_parameters['resolution']}

    def get_rf_cavity_system(self, twiss=None):
        '''
        Obtain RF-cavity related quantities like energy acceptance and synchrotron tune.

        The underlying MAD-X model is given by
        V_RF = V*sin(2*pi(lag - harmonic*f0*t)) ,
        where V is given in terms of MV.
        '''
        if twiss == None:
            if self.verbose:
                print ('get_rf_cavity_system: Performing MAD-X twiss on given lattice ...')
            twiss = self.twiss()
        indices_rfcavities = twiss.keyword == 'rfcavity'
        voltages = twiss.volt[indices_rfcavities]*1e6
        harmonics = twiss.harmon[indices_rfcavities]
        phase_lags = twiss.lag[indices_rfcavities]*2*np.pi
        frequencies = twiss.freq[indices_rfcavities]
        return {'voltage': voltages, 'harmonic': harmonics, 'phase_rad': phase_lags, 'frequency': frequencies}

    def energy_acceptance(self, phase, harmonic, voltage, slip_factor):
        '''
        Compute energy acceptance and synchrotron revolution of a single RF cavity.
        '''
        # compute energy acceptance (see Frank Tecker's notes https://arxiv.org/pdf/2011.02932.pdf)
        G_phi = 2*np.cos(phase) + (2*phase - np.pi)*np.sin(phase)
        beta0 = self.beam.beta.value
        charge = constants.e
        energy0_SI = self.beam.energy.value*1e9*constants.e
        return beta0*np.sqrt(-charge*voltage/(np.pi*harmonic*slip_factor*energy0_SI)*G_phi)

    def update_beam(self, natural_parameters: dict, update_ey=False, ey_scale=1):
        '''
        Update beam parameters using the results stored in the dictionary natural parameters.

        update_ey: If True, then also set the y-emittance according to the theoretical lower limit.
        ey_scale: (default 1) scaling factor for the lower bound of the natural y-emittance. Only in effect
        if update_ey == True.
        '''
        old_ex_value = self.beam.ex.value
        self.beam.ex.value = natural_parameters['nat_ex']
        old_exn_value = self.beam.exn.value
        self.beam.exn.value = self.beam.ex.value*self.beam.gamma.value

        if update_ey:
            old_ey_value = self.beam.ey.value
            self.beam.ey.value = natural_parameters['nat_ey_lb']*ey_scale
            old_eyn_value = self.beam_eyn.value
            self.beam.eyn.value = self.beam.ey.value*self.beam.gamma.value

        old_sigt_value = self.beam.sigt.value
        self.beam.sigt.value = natural_parameters['nat_z']
        old_sige_value = self.beam.sige.value
        self.beam.sige.value = natural_parameters['nat_dee']

        if self.verbose:
            print ('\nUpdating beam parameters using natural values')
            print (f'ex: OLD: {old_ex_value}, NEW: {self.beam.ex.value}')
            print (f'exn: OLD: {old_exn_value}, NEW: {self.beam.exn.value}')
            if update_ey:
                print (f'ey: OLD: {old_ey_value}, NEW: {self.beam.ey.value}')
                print (f'eyn: OLD: {old_eyn_value}, NEW: {self.beam.eyn.value}')
            print (f'sigt: OLD: {old_sigt_value}, NEW: {self.beam.sigt.value}')
            print (f'sige: OLD: {old_sige_value}, NEW: {self.beam.sige.value}')

    def touschek_lifetime(self, precise=False, precision=16, symmetry=1, **kwargs):
        '''
        Compute the touschek lifetime for the given optics.

        Input parameters:
        see touschek.touschek.lifetime
        '''
        nat = self.get_natural_parameters(**kwargs)
        self.update_beam(nat)
        return lifetime(optics=self, delta_pm=nat['momentum_acceptance'], 
                        precise=precise, precision=precision, symmetry=symmetry, verbose=self.verbose)

    def plot_survey(self, **kwargs):
        return plot_survey(madx=self.madx, **kwargs)

    def plot_touschek_losses(self, touschek_results, **kwargs):
        plot_touschek_losses(optics=self, touschek_results=touschek_results, **kwargs)

    def test_tune(self, **kwargs):
        if not hasattr(self, 'function'):
            raise ValueError('Error: run compute_optics_functions first.')
        test_tune(optics=self, **kwargs)

    def track(self, coordinates, turns: int, **kwargs):
        '''
        Track a set of particles through the optics.
        '''
        self.madx_thin.track(**kwargs)
        coordinate_keys = ['x', 'px', 'y', 'py', 't', 'pt']
        n_particles = len(coordinates)
        for point in coordinates:
            self.madx_thin.start(x=point[0], px=point[1], y=point[2], py=point[3], t=point[4], pt=point[5])
        self.madx_thin.run(turns=turns)
        self.madx_thin.endtrack()
        coordinates = np.array([self.madx_thin.table.tracksumm[ck][n_particles:] for ck in coordinate_keys]).T
        return coordinates

    def track_series(self, coordinates, turns, **kwargs):
        '''
        Track a series of turns; the output coordinates for each series are taken as input arguments for the next.
        '''
        coordinates_all = []
        for t in np.diff([0] + list(turns)):
            coordinates = self.track(coordinates=coordinates, turns=t, **kwargs)
            coordinates_all.append(coordinates)
        return np.vstack(coordinates_all)
