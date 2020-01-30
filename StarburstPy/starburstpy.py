# -*- coding: utf-8 -*-
"""
    Main code for StarburstPy.
    
    See example scripts for how to use.

    author: Ryan Tanner
    email: ryan.tanner@nasa.gov

"""
import StarburstPy as sb
import os
import numpy as np
import subprocess

class SbInput(object):
    
    def __init__(self, model_name, chatter = 2):
        """
        Creates starburst input object.
        
            - A model name is required.
        
            - Optional: chatter, this controls how much output is printed to stdout.
                chatter = 0 : No output
                chatter = 1 : Only errors
                chatter = 2 : Errors and warnings [default]
                chatter = 3 : Errors, warnings, and messages
        """
        
        self.model_name = model_name
        self.sb_messages = sb.sb_messages(chatter = chatter)
        self.init_default()
        
        
    def init_default(self):
        """
        Sets default values for input when object is created.
        """
        self.hi_res_set = False
        self.set_star_formation()
        self.set_IMF()
        self.set_SN_cut_off()
        self.set_evolve_track()
        self.set_wind_model()
        self.set_time()
        self.set_mass_grid()
        self.set_mass_track()
        self.dt_print_spectra()
        self.set_atmo()
        self.high_res_Z()
        self.uv_line_lib()
        self.rgs_turb_abund()
        self.set_output()
        
#        self.sb_messages.message('Total mass set to default {0} M_sun'.format(self.total_mass),'SbInput.init_default')
        """
        TO DO: Print out all defaults as messages.
        """
        
        
    def set_model_name(self, model_name = None):
        """
        The name for creating the output files.
        Usually the model name is passed in as an argument when SbInput is called,
        but this function can be used to create a different set of output files
        if running multiple models.
        """
        
        if model_name is None:
            self.model_name = 'Default_Model'
        else:
            self.model_name = model_name
        
        
    def set_star_formation(self, fixed_mass = True, total_mass = 1e6, SFR = 1.0):
        """
        Function for setting either the star formation rate for continuous SF, or
        using an instantaneous burst with a total fixed mass.
        
            fixed_mass: [boolean] True sets an instantaneous burst with a total fixed mass.
                                  False sets continuous SF with a constant SFR.
            total_mass: [number]  In solar masses, the total mass of an inst. burst 
            SFR:        [number]  In solar masses/year, star formation rate.
        """
        
        self.fixed_mass = fixed_mass
        self.total_mass = total_mass
        self.SFR = SFR
        if fixed_mass:
            self.sb_messages.message('Total mass set to {0} M_sun'.format(self.total_mass),'SbInput.set_star_formation')
        else:
            self.sb_messages.message('SFR set to {0} M_sun/yr'.format(self.SFR),'SbInput.set_star_formation')
                            
    
    def set_IMF(self, IMF = None, IMF_exp = None, IMF_bound = None, mass_lo = None, mass_hi = None):
        """
        Function to set the IMF. A standard IMF can be chosen, or a specialized one can be entered.
            User can choose a preset: IMF = 'Kroupa' or IMF = 'Salpeter'
            Or user can set their own IMF exponents and bounds
            
            IMF: [str] Chose a preset IMF. 
                        -- Kroupa [default]
                        -- Salpeter
            IMF_exp: [list] List of IMF exponents, e.g. [1.3,2.3]
            IMF_bound: [list] IMF mass boundaries. The number of boundaries MUST be one more
                       than the number of exponents. Ex. if IMF_exp = [1.3,2.3] then
                       IMF_bound could be [0.1, 0.5, 120.0], where the first (last) number
                       is the lowest (highest) mass limit of stars in the model. The middle
                       number is the dividing point between the two exponents.
            mass_lo: [number] Lower mass limit for default IMF's. Default is 0.08 M_sun.
            mass_hi: [number] Upper mass limit for default IMF's. Default is 120.0 M_sun.
        """
        
        if mass_lo is None:
            self.mass_lo = 0.08
            
        if mass_hi is None:
            self.mass_hi = 120.0
        
        if IMF_exp is not None and IMF_bound is not None:
            self.IMF_exp = IMF_exp
            self.IMF_bound = IMF_bound
        else:
            if IMF == 'Kroupa':
                self.IMF = IMF
                self.IMF_exp = [1.3, 2.3]
                self.IMF_bound = [self.mass_lo, 0.5, self.mass_hi]
            elif IMF == 'Salpeter':
                self.IMF = IMF
                self.IMF_exp = [2.35]
                self.IMF_bound = [self.mass_lo, self.mass_hi]
            elif IMF is None:
                self.IMF = 'Kroupa'
                self.IMF_exp = [1.3, 2.3]
                self.IMF_bound = [self.mass_lo, 0.5, self.mass_hi]
                self.sb_messages.message('IMF set to default Kroupa.','SbInput.set_IMF')
            else:
                self.sb_messages.error('IMF not recognized.','SbInput.set_IMF')
                
        self.IMF_exp = np.array(self.IMF_exp)
        self.IMF_bound = np.array(self.IMF_bound)
                
                
    def set_SN_cut_off(self, lo_mass = 8.0, hi_mass = 120.0):
        """
        Function to set supernova mass limits in M_sun.
            lo_mass: [number] Lower mass limit, default 8.0 M_sun
            hi_mass: [number] Upper mass limit, default 120 M_sun
        """
        
        self.SN_lo_mass = lo_mass
        self.SN_hi_mass = hi_mass
            
            
    def set_evolve_track(self, track = None, Z = None):
        """
        Function for setting the evolutionary track and metallicity.
            track: [str] Sets evolutionary track
                Options:
                    PadovaOrig
                    PadovaAGB
                    GenevaStd
                    GenevaHigh
                    Genevav00 [default] Z = 0.014 (Solar)
                    Genevav40
            Z: [number] Metallicity
                Your choice of evolutionary track will determine which
                metallicites you can use.
                    Padova: [0.0004, 0.004, 0.008, 0.02, 0.05]
                    Geneva(Std & High): [0.001, 0.004, 0.008, 0.02, 0.04]
                    Geneva(v00 & v40): [0.001, 0.002, 0.008, 0.014, 0.04]
                        Z = 0.001, 0.008, and 0.04 are not available.
        """
        
        tracks = ['PadovaOrig', 'PadovaAGB', 'GenevaStd', 'GenevaHigh', 'Genevav00', 'Genevav40']
        
        if track is None:
            self.evolve_track = 'Genevav00'
            self.sb_messages.message('Evolutionary track set to Genevav00.','SbInput.set_evolve_track')
        elif any(x == track for x in tracks):
            self.evolve_track = track
        else:
            self.sb_messages.error('Evolutionary track not recognized. Available options {0}'.format(tracks),'SbInput.set_evolve_track')

        PZ = [0.0004, 0.004, 0.008, 0.02, 0.05]
        GSH = [0.001, 0.004, 0.008, 0.02, 0.04]
        Gv = [00000, 0.002, 00000, 0.014, 0000]
        
        """
        Corresponding mass limits for WR star formation.
        """
        PZ_mass_lim = [61, 42, 35, 25, 21]
        GS_mass_lim = [80, 52, 42, 32, 25]
        GH_mass_lim = [61, 42, 35, 25, 21]
        v0_mass_lim = [84, 84, 25, 25, 25]
        v4_mass_lim = [55, 55, 20, 20, 20]

        if Z is None:
            self.metallicity = 0.014
            self.sb_messages.message('Metallicity set to default {0}.'.format(self.metallicity),'SbInput.set_evolve_track')
        else:
            if self.evolve_track == 'PadovaOrig' or self.evolve_track == 'PadovaAGB':
                if any(x == Z for x in PZ):
                    self.metallicity = Z
                    ind = PZ.index(Z)
                    self.wr_mass_lim = PZ_mass_lim[ind]
                else:
                    self.sb_messages.error('Metallicity not allowed. Available options for {0} are {1}'.format(self.evolve_track,PZ),'SbInput.set_evolve_track')
            elif self.evolve_track == 'GenevaStd' or self.evolve_track == 'GenevaHigh':
                if any(x == Z for x in GSH):
                    self.metallicity = Z
                    ind = GSH.index(Z)
                    if self.evolve_track == 'GenevaStd':
                        self.wr_mass_lim = GS_mass_lim[ind]
                    if self.evolve_track == 'GenevaHigh':
                        self.wr_mass_lim = GH_mass_lim[ind]
                else:
                    self.sb_messages.error('Metallicity not allowed. Available options for {0} are {1}'.format(self.evolve_track,GSH),'SbInput.set_evolve_track')
            elif self.evolve_track == 'Genevav00' or self.evolve_track == 'Genevav40':
                if any(x == Z for x in Gv):
                    self.metallicity = Z
                    ind = Gv.index(Z)
                    if self.evolve_track == 'Genevav00':
                        self.wr_mass_lim = v0_mass_lim[ind]
                    if self.evolve_track == 'Genevav40':
                        self.wr_mass_lim = v4_mass_lim[ind]
                else:
                    self.sb_messages.error('Metallicity not allowed. Available options for {0} are {1}'.format(self.evolve_track,Gv),'SbInput.set_evolve_track')

        if self.hi_res_set:
            self.high_res_Z(Z = self.metallicity)
            

    def set_wind_model(self, wind = None):
        """
        Selects model to calculate wind power.
            wind: [str] Possible values 
                    -- Evolution[default]
                    -- Emperical
                    -- Theoretical
                    -- Elson
        """
        
        models = ['Evolution', 'Emperical', 'Theoretical', 'Elson']
        
        if wind is None:
            self.wind_model = 'Evolution'
        elif any(x == wind for x in models):
            self.wind_model = wind
        else:
            self.sb_messages.error('Wind model not recognized. Available options are {0}'.format(models),'SbInput.set_wind_model')
        
        
    def set_time(self, start_time = 0, end_time = 5e7, log_scale = False, nsteps = 1000, dt = 1e5):
        """
        Sets parameters for the time integration.
            -- start_time: Initial time in yrs. Note SB99 warned that start time 
                           must be small but not zero. If set to zero then it will be 
                           automatically set to a small number.
            -- end_time:   Final time in yrs.
            
            -- log_scale:  [boolean] If True, then integration uses timesteps set
                           logarithmically. If False, then linear timesteps used.
            -- nsteps:     [number] Number of time steps to take.
                           Only used if log_scale = True
            -- dt:         [number] Size of time step in yrs.
                           Only used if log_scale = False
        """
        
        if start_time <= 0:
            start_time = 1e4
            self.sb_messages.message('Start time set to {0} yrs.'.format(start_time),'SbInput.set_time')
        
        self.start_time = start_time
        self.end_time = end_time

        #if log_scale is true then time steps will be taken logarithmically
        
        self.log_scale = log_scale
        
        self.dt = dt # Only used if log_scale = False
        self.nsteps = nsteps # Only used if log_scale = True


    def set_mass_grid(self, mass_grid = None):
        """
        Set mass grid.
            mass_grid: [str] Possible values are:
                        -- Full_Isochrone [default]
                        -- Isochrone_Large
                        -- Large
                        -- Small
        """
        
        grids = ['Full_Isochrone', 'Isochrone_Large', 'Large', 'Small']
        
        if mass_grid is None:
            self.mass_grid = 'Full_Isochrone'
        elif any(x == mass_grid for x in grids):
            self.mass_grid = mass_grid
        else:
            self.sb_messages.error('Mass grid not recognized. Available options are {0}'.format(grids),'SbInput.set_mass_grid')
        
        
    def set_mass_track(self, mass_range = None):
        """
        For debugging and special use. See original note in Starburst99.
        For now set to null [0,0] by default. RT 12/17/2019
        """
        
        if mass_range is None:
            self.mass_track = [0,0]
        else:
            if len(mass_range) == 2 and mass_range[0] <= mass_range[1]:
                self.mass_track = mass_range
            else:
                self.sb_messages.error('Mass range not allowed.','SbInput.mass_track')
    
    
    def dt_print_spectra(self, dt = 1e6):
        """
            dt: [number] Time step for printing out spectra.
        """

        self.spec_dt = dt
        
            
    def set_atmo(self, atmo = None):
        """
        Selects model atmosphere for low resolution spectra output.
            atmo: [str] Possible values:
                        -- Planck
                        -- Lejeune
                        -- Lejuene-Schmutz
                        -- Lejeune-Hiller
                        -- Pauldrach-Hiller [default]
        """
        
        models = ['Planck', 'Lejeune', 'Lejuene-Schmutz', 'Lejeune-Hiller', 'Pauldrach-Hiller']
        
        if atmo is None:
            self.low_res_spec = 'Pauldrach-Hiller'
        elif any(x == atmo for x in models):
            self.low_res_spec = atmo
        else:
            self.sb_messages.error('Atmosphere model not recognized. Available options are {0}'.format(models),'SbInput.set_atmo')
        

    def high_res_Z(self, Z = None):
        """
        Sets metallicity for the optical high res spectra library.
            Options are: Z = 
                            -- 0.001
                            -- 0.008
                            -- 0.02 [default]
                            -- 0.04
        """
        
        warn = False
        
        if Z == 0.0004:
            Z = 0.001
            warn = True
        elif Z == 0.004 or Z == 0.002:
            Z = 0.001
            warn = True
        elif Z == 0.014:
            Z = 0.02
            warn = True
        elif Z == 0.050:
            Z = 0.04
            warn = True
        
        if warn:
            self.sb_messages.warning('Metallicity for optical high res library set to {0}'.format(Z),
                                     'SbInput.high_res_Z')
        
        Zs = [0.001, 0.008, 0.02, 0.04]
        
        if Z is None:
            self.high_res_spec_Z = Zs[2]
        elif any(x == Z for x in Zs):
            self.high_res_spec_Z = Z
        else:
            self.sb_messages.error('Metallicity not allowed. Available options for high res spectra are {1}'.format(Zs),
                                   'SbInput.high_res_Z')
        
        self.hi_res_set = True
        
        
    def uv_line_lib(self, lib = None):
        """
            UV_line_lib: [str] Choice of the UV spectral library:
                            -- Solar [default]
                            -- LMC/SMC
        """
        
        libs = ['Solar', 'LMC/SMC']
        
        if lib is None:
            self.UV_line_lib = libs[0]
        elif any(x == lib for x in libs):
            self.UV_line_lib = lib
        else:
            self.sb_messages.error('UV line library not allowed. Available options are {1}'.format(libs),'SbInput.UV_line_lib')


    def rgs_turb_abund(self, RGS_turb_vel = 3, RGS_abund = True):
        """
        RGS microturbulance and abundance.
            
            Possible values:
                0 <= RGS_turb_vel <= 6
                RGS_abund = True/False
        """

        if RGS_turb_vel >= 0 and RGS_turb_vel <= 6:
            self.RGS_turb_vel = RGS_turb_vel
        else:
            self.sb_messages.error('Microturbulance value not allowed. Must be 0 <= vel <= 6.','SbInput.rgs_turb_abund')
            
        self.RGS_abund = RGS_abund
        
        
    def set_output(self, writeinput = False, silent = False, QUANTA = True, SNR = True, HRD = False, POWER = True, 
                   SP = True, YIELDS = True, SPECTRUM = True, LINE = True, 
                   COLOR = True, WIDTH = True, FEATURES = True, OVI = True, 
                   HIRES = True, WRLINES = True, IFASPEC = True,
                   out_type = None):
        """
            -- writeinput -- Controls whether input parameters are written to a
                             textfile, separate from the model_name.input file 
                             written for Starburst99.
                                
            -- silent     -- Surpresses all file output. Data only kept in out_data 
                             object that can be used in python. Use this if you
                             only want to work with the output data in an active 
                             python session.
                             Overrides all other output choices. 
                             
            -- out_type   -- Output filetypes.
                Possible values:
                    - original -- [default] Keeps files in the original Starburst99 file format.
                    - hdf5 -- Converts the original Starburst99 files into hdf5 file format.
                    
        Output files controled by True/False. By default all are True, except HRD (3).
        
        OUTPUT FILES:
            
            0 model_name.input    A listing of all input parameters.
            1 model_name.QUANTA   Computation of the number of ionizing photons. Depends on 7.
            2 model_name.SNR      Calculation of the supernova rate and the mechanical luminosities (winds and supernovae). It requires (4) to obtain the stellar wind luminosities.
            3 model_name.HRD      HRD with a few evolutionary tracks. Off by default.
            4 model_name.POWER    Mechanical luminosity and related quantities due to winds.
            5 model_name.SP       Two output files containing the stellar spectral types during each time step and the relative numbers of WR stars.
            6 model_name.YIELDS   The mass in individual elements released via stellar winds and supernovae.
            7 model_name.SPECTRUM The continuous spectrum of the stellar population for each time step.
            8 model_name.LINE     The ultraviolet line spectrum at 0.75 A resolution from 1200 to 1600 A (LMC/SMC library) or to 1800 A (Milky Way library). Depends on 7 and 1.
            9 model_name.COLOR    Calculation of colors and magnitudes. Depends on 7 and 1.
           10 model_name.WIDTH    Calculation of the strengths of H_alpha, H_beta, Pa_beta, and Br_gamma. Depends on 7 and 1.
           11 model_name.FEATURES Calculation of the strengths of various spectral features. 
           12 model_name.OVI      This subroutine is equivalent to (8), but it computes the spectral region between 1000 and 1180 A.
           13 model_name.HIRES    Calculates fully theoretical spectra between 3000 and 7000 A at 0.3 A resolution.
           14 model_name.WRLINES  Calculation of the most important WR emission lines. Depends on 7 and 1.
           15 model_name.IFASPEC  Calculation of a high-resolution UV line spectrum from model atmospheres, as opposed to using an empirical library.
        
        Dependencies: (xxx --> depends on yyy)
        
            1 --->  7
            2 --->  4
            8 --->  1 --->  7
            9 --->  1 --->  7
           10 --->  1 --->  7
           14 --->  1 --->  7
           
           Note: 8, 9, 10, and 14 can be output without 1 but not without 7, as long as 
                 the effect from 1 is small.
        
        """
        
        """
        Check dependencies.
        """
        
        if QUANTA and not SPECTRUM:
            self.sb_messages.warning('Setting output SPECTRUM = True. \n Computation of the number of ionizing photons depends on the continuous spectrum.','set_output')
            SPECTRUM = True
        if SNR and not POWER:
            self.sb_messages.warning('Setting output POWER = True. \n Calculation of the mechanical luminosities depends related quantities due to winds.','set_output')
            POWER = True
        if LINE and not SPECTRUM:
            self.sb_messages.warning('Setting output SPECTRUM = True. \n The ultraviolet line spectrum depends on the continuous spectrum.','set_output')
            SPECTRUM = True
        if COLOR and not SPECTRUM:
            self.sb_messages.warning('Setting output SPECTRUM = True. \n Calculation of colors and magnitudes depends on the continuous spectrum.','set_output')
            SPECTRUM = True
        if WIDTH and not SPECTRUM:
            self.sb_messages.warning('Setting output SPECTRUM = True. \n Calculation of line strengths depends on the continuous spectrum.','set_output')
            SPECTRUM = True
        if WRLINES and not SPECTRUM:
            self.sb_messages.warning('Setting output SPECTRUM = True. \n Calculation of the most important WR emission lines depends on the continuous spectrum.','set_output')
            SPECTRUM = True
        
        if LINE and not QUANTA:
            self.sb_messages.warning('It is advised to set QUANTA = True if LINE = True.','set_output')
        if COLOR and not QUANTA:
            self.sb_messages.warning('It is advised to set QUANTA = True if COLOR = True.','set_output')
        if WIDTH and not QUANTA:
            self.sb_messages.warning('It is advised to set QUANTA = True if WIDTH = True.','set_output')
        if WRLINES and not QUANTA:
            self.sb_messages.warning('It is advised to set QUANTA = True if WRLINES = True.','set_output')
        
        self.outputs = {}
        
        self.outputs['input'] = writeinput
        self.outputs['QUANTA'] = QUANTA
        self.outputs['SNR'] = SNR
        self.outputs['HRD'] = HRD
        self.outputs['POWER'] = POWER
        self.outputs['SP'] = SP
        self.outputs['YIELDS'] = YIELDS
        self.outputs['SPECTRUM'] = SPECTRUM
        self.outputs['LINE'] = LINE
        self.outputs['COLOR'] = COLOR
        self.outputs['WIDTH'] = WIDTH
        self.outputs['FEATURES'] = FEATURES
        self.outputs['OVI'] = OVI
        self.outputs['HIRES'] = HIRES
        self.outputs['WRLINES'] = WRLINES
        self.outputs['IFASPEC'] = IFASPEC
        
        if out_type is None:
            self.outputs['out_type'] = 'original'
        elif out_type == 'hdf5':
            self.outputs['out_type'] = out_type
        elif out_type == 'original':
            self.outputs['out_type'] = out_type
        else:
            self.sb_messages.error('File type for output files not recognized. Possible types are original or hdf5.','set_output')
        
        if silent:
            self.outputs['input'] = False
            self.outputs['QUANTA'] = False
            self.outputs['SNR'] = False
            self.outputs['HRD'] = False
            self.outputs['POWER'] = False
            self.outputs['SP'] = False
            self.outputs['YIELDS'] = False
            self.outputs['SPECTRUM'] = False
            self.outputs['LINE'] = False
            self.outputs['COLOR'] = False
            self.outputs['WIDTH'] = False
            self.outputs['FEATURES'] = False
            self.outputs['OVI'] = False
            self.outputs['HIRES'] = False
            self.outputs['WRLINES'] = False
            self.outputs['IFASPEC'] = False
        
        
    def run_starburst(self, run_SB99 = True):
        """
        Method to run the starburst model. Simply calls the function "run_starburst" defined below.
        """
        output_data = sb.out_data(self.model_name)
        output_data = run_starburst(self, run_SB99 = run_SB99)
        return output_data
        
        
class out_data(object):
    """
        Class for holding all output data. Contains a method for reading in data
        from previous runs.
        
        Contains four dictionaries:
            -- data        : Contains all output data in sub-dictionaries.
                             Each entry in 'data' is a dictionary for each output file.
                             Columns of data are stored as numpy arrays.
                             
            -- headers     : Contains headers from output files.
            
            -- fname2dname : Contains the original file endings and 
                             corresponding dictionary key in 'data'.
                             
            -- data_dict   : Each item in the dictionary corresponds to an output
                             file and contains a list of names of numpy arrays in
                             the 'data' dictionary.
    """
    
    def __init__(self, model_name, output_dir = None):
        
        self.model_name = model_name
        if output_dir is None:
            self.output_dir = sb.indata.output_dir
        else:
            self.output_dir = output_dir
        self.data = {}
        self.headers = {}
        self.fname2dname = {'color' : 'color',
                    'ewidth' : 'recombination_lines',
                    'hires' : 'hires_spectra',
                    'ifaspec' : 'infrared_spectra',
                    'ovi' : 'uv_spectra1',
                    'power' : 'power',
                    'quanta' : 'quanta',
                    'snr' : 'snr',
                    'spectrum' : 'spectrum',
                    'sptyp1' : 'spectral_type1',
                    'uvline' : 'uv_spectra2',
                    'yield' : 'chemical_yields'}
        self.data_dict = {'color' : ['Time', '130-V', '210-V', 'U-B', 'B-V', 'V-R', 'V-I', 
                               'V-J', 'V-H', 'V-K', 'V-L', 'Mag_V', 'Mag_B', 'Mag_Bol'],
                    'recombination_lines' : ['Time', 'Continuum_H_A', 'Luminosity_H_A', 'Equ_width_H_A', 
                                'Continuum_H_B', 'Luminosity_H_B', 'Equ_width_H_B', 
                                'Continuum_P_B', 'Luminosity_P_B', 'Equ_width_P_B', 
                                'Continuum_B_G', 'Luminosity_B_G', 'Equ_width_B_G'],
                    'hires_spectra' : ['Time', 'Wavelength', 'Log_Luminosity', 'Normalized'],
                    'infrared_spectra' : ['Time', 'Wavelength', 'Log_Luminosity', 'Normalized'],
                    'uv_spectra1' : ['Time', 'Wavelength', 'Log_Luminosity', 'Normalized'],
                    'power' : ['Time', 'Power_All', 'Power_OB', 'Power_RSG','Power_LBV','Power_WR',
                               'Energy','Momentum_All', 'Momentum_OB', 'Momentum_RSG','Momentum_LBV','Momentum_WR'],
                    'quanta' : ['Time', 'HI_rate', 'HI%', 
                                'HeI_rate', 'HeI%',
                                'HeII_rate', 'HeII%', 'Log_Luminosity'],
                    'snr' : ['Time', 'All_Rate', 'All_Power', 'All_Energy', 
                             'Type_IB_Rate', 'Type_IB_Power', 'Type_IB_Energy', 
                             'Typical_Mass', 'Lowest_Mass', 
                             'Stars_SN_Power', 'Stars_SN_Energy'],
                    'spectrum' : ['Time', 'Wavelength', 'Total', 'Stellar','Nebular'],
                    'spectral_type1' : ['Time', 'All','O','WR','WN','WC','WR/O','WC/WN'],
                    'uv_spectra2' : ['Time', 'Wavelength', 'Log_Luminosity', 'Normalized'],
                    'chemical_yields' : ['Time','H','He','C','N','O','Mg','Si','S','Fe',
                               'All_winds','All_SN','All_winds_SN','Total_mass']}
        
     
    def read_data_files(self, model_name = None, output_dir = None):
        """
        Method for reading in Starburst99 output. This can be used to load existing 
        data.
        """
        
        if model_name is None:
            model_name = self.model_name
            
        if output_dir is None:
            output_dir = self.output_dir
        
        self.data, self.headers = sb.read_output_data(model_name = model_name, output_dir = output_dir)
        
        
def run_starburst(input_param, run_SB99 = True):
    """
    This is where the work happens. SbInput object must be passed in.
    
    Organized as follows:
        
        1. Write input parameters to file.
        
        2. Run Starburst99.
        
        3. Read output files written by Starburst99. 
           Convert to hdf5 if necessary. [Not implemented yet.]
           
    """
    output_data = sb.out_data(input_param.model_name)
    
    """
    Step 1
    """
    
    output_dir = sb.indata.output_dir
    file_name = input_param.model_name+'_input.txt'
    full_name = os.path.join(output_dir,file_name)
    if os.path.exists(full_name):
        os.remove(full_name)
    if input_param.outputs['input']:
        sb.write_input_txt(input_param, full_name)
    
    """
    Step 2
    """
    
    if run_SB99:
        file_name = input_param.model_name+'.input'
        full_name = os.path.join(output_dir,file_name)
        sb.write_input_original(input_param, full_name)
        output_data.data, output_data.headers = run_original_SB99(input_param.model_name)
        
    """
    Step 3
    """

    return output_data


def run_original_SB99(model_name):
    """
    Runs the original Fortran77 code of Starburst99. Creates a 'go_galaxy.sh' script
    to run the code. This is the same script that comes with the code.
    
    Output directory and Starburst99 directory need to be set previously in 
    sb.indata.output_dir and sb.indata.SB99_dir respectively.
    """
    
    output_dir = sb.indata.output_dir
    rel_path = os.path.relpath(sb.indata.SB99_dir,output_dir)
    file_name = model_name+'.input'
    full_name = os.path.join(output_dir,file_name)
    if not os.path.exists(full_name):
        sb.sb_mess.error('Input File does not exist in {0}.'.format(output_dir),'run_original_SB99')
    run_file = os.path.join(output_dir,'go_galaxy.sh')
    if os.path.exists(run_file):
        os.remove(run_file)
    if not os.path.exists(os.path.join(sb.indata.SB99_dir,'galaxy')):
        sb.sb_mess.error('Starburst99 exacutable not found in {0} directory.'.format(sb.indata.SB99_dir),'run_original_SB99')
        
    f = open(run_file, 'w')
    
    lines = ['#! /bin/csh -f\n',
            '# -----------------------------------------\n',
            '# UNIX Script to run galaxy code\n',
            '# -----------------------------------------\n',
            '#\n',
            '# Usage: 1) Edit six user-defined quantities found below\n',
            '#        2) Run code with \'nice go_galaxy &\' --> output in \'time_used\'\n',
            '# -----------------------------------------------------------------------\n',
            '# Six user-defined quantities:\n',
            '#           Directory where the run is made (=directory of this file):\n',
            'set    drun={0}\n'.format(output_dir),
            '#           Name of input file:\n',
            'set  ninput={0}\n'.format(model_name+'.input'),
            '#           Name and extension number of output files:\n',
            '#           --> files will be: noutput.colornext, noutput.quantnext etc.\n',
            'set noutput={0}\n'.format(model_name),
            'set    next=1\n',
            '#           Directory where code is:\n',
            'set   dcode={0}\n'.format(rel_path),
            '#          Name of code: \n',
            'set   ncode=galaxy\n',
            '#           Directory where libraries are:\n',
            'set   dlib={0}\n'.format(rel_path),
            '#- - - - - - - DO NOT MODIFY below this line - - - - - - - - - - - - - -\n',
            '\n',
            'cd $drun\n',
            '# A) Create links (assign) to directories with input files\n',
            '	# Tracks + Spectral type calibration:\n',
            'if (!(-e tracks)) ln -s $dlib/tracks/ tracks\n',
            '	# Atmosphere models:\n',
            'if (!(-e lejeune)) ln -s $dlib/lejeune/ lejeune\n',
            '	# Spectral libraries:\n',
            'if (!(-e auxil)) ln -s $dlib/auxil/ auxil\n',
            '\n',
            '# B) Link input file, RUN code, and save output data\n',
            'echo "Job on `hostname` started at `date`" >time_used\n',
            'echo " ">>time_used\n',
            '\n',
            'if (-e fort.1) /bin/rm fort.1\n',
            'if (-e $ninput) ln -s $ninput fort.1\n',
            '(/usr/bin/time $dcode/$ncode) >>& time_used\n',
            'echo " ">>time_used\n',
            ' \n',
            '$dcode/save_output $noutput $next >>time_used\n',
            'echo "Done at `date`" >>time_used\n']
    
    f.writelines(lines)
    
    f.close()
    
    subprocess.call(['/bin/csh','-c',run_file])
    
    data, headers = sb.read_output_data(model_name = model_name, output_dir = output_dir)
    
    return data, headers
    