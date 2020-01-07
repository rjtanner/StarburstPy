# -*- coding: utf-8 -*-

import StarburstPy as sb

"""

    A script that shows all the options for StarburstPy.

*********************** TOO MUCH? ***********************

    Try run_sbpy_simp.py. Only the most common parameters 
    are shown in that example script.

    Author: Ryan Tanner
    email: ryan.tanner@nasa.gov

"""

# Sets path links to current directory, output directory, and the directory where
# the Starburst99 fortran code is compiled. Defaults set to current directory.

datapaths = sb.indata
# datapaths.cur_dir is available for the current directory.
datapaths.output_dir = './sboutput' # Directory where output files should go.
                                    # If the directory does not exist it will be
                                    # created.
datapaths.SB99_dir = '../../starburst99' # Either absolute or relative path to 
                                         # the directory where Starburst99 is 
                                         # compiled.

# A model name must be passed in to make the input object.

model_name = 'happy'

# The input object "sbinput" contains all the input parameters and methods for 
# setting all parameters, and a method for running Starburst99.
# If just the model name is given then the default parameters are set to the
# defaults from the Starburst99 webpage.
#
# rank sets how much is written to standard output.
# rank = 0 : No output
# rank = 1 : Only errors
# rank = 2 : Errors and warnings [default]
# rank = 3 : Errors, warnings, and messages

sbinput = sb.SbInput(model_name, rank = 3)

# fixed_mass = True/False
#   If True: Star formation with a fixed mass.
#   If False: Continuous star formation.
#
# Based on your choice of fixed_mass either the total_mass or the SFR needs
# to be set.
# total_mass in M_sun. Only used if fixed_mass = True
# SFR in M_sun/yr. Only used if fixed_mass = False

sbinput.set_star_formation(fixed_mass = False, SFR = 1.0)
# The next line will change the settings from the previous line.
sbinput.set_star_formation(fixed_mass = True, total_mass = 1e6)

# Preset IMF available 'Kroupa', 'Salpeter'
# IMF_exp: IMF exponents [a1, a2, a3, ....]
# IMF_bound: IMF mass boundaries over which the exponents [M1, M2, M3, ...]
#   where M1->M2 is the mass range for the exponent a1.
#   The number of mass boundaries needs to be the number of exponents +1
#
# Note: The number of IMF intervals does not need to be entered.
# 
# mass_lo: Lowest mass in the IMF (note some later choices may make these settings moot).
# mass_hi: Highest mass in the IMF.

sbinput.set_IMF(IMF = 'Salpeter')
# The next line will change the settings from the previous line.
sbinput.set_IMF(IMF = 'Kroupa')
# The next line will change the settings from the previous line.
sbinput.set_IMF(IMF_exp = [1.3,2.3], IMF_bound = [0.1, 0.5, 120.0], mass_lo = 0.1, mass_hi = 120.0)

# lo_mass: Lowest mass that stars will go supernova.
# hi_mass: Highest mass that stars will go supernova. Form black holes instead 
#          and do not contribute to mass and energy feedback. 

sbinput.set_SN_cut_off(lo_mass = 8.0, hi_mass = 120.0)

# track: Sets evolutionary track
#          Options:
#            PadovaOrig
#            PadovaAGB
#            GenevaStd
#            GenevaHigh
#            Genevav00 [default] Z = 0.014 (Solar)
#            Genevav40
# Z: Metallicity
#          Your choice of evolutionary track will determine which
#          metallicites you can use.
#            Padova: [0.0004, 0.004, 0.008, 0.02, 0.05]
#            Geneva(Std & High): [0.001, 0.004, 0.008, 0.02, 0.04]
#            Geneva(v00 & v40): [0.001, 0.002, 0.008, 0.014, 0.04]
#              NOTE: Z = 0.001, 0.008, and 0.04 are not available for Geneva(v00 & v40)

sbinput.set_evolve_track(track = 'Genevav40', Z = 0.002)

#Selects model to calculate wind power.
#    wind: Possible values [Evolution, Emperical, Theoretical, Elson]

sbinput.set_wind_model(wind = 'Evolution')

# Sets parameters for the time integration.
#   start_time: Initial time in yrs. Note SB99 warned that start time 
#               must be small but not zero. If set to zero then it will be 
#               automatically set to a small number.
#
#   end_time:   Final time in yrs.
#            
#   log_scale: If True, then integration uses timesteps set
#              logarithmically. If False, then linear timesteps used.
#
#   nsteps:    Number of time steps to take.
#              Only used if log_scale = True
#
#   dt:        Size of time step in yrs.
#              Only used if log_scale = False

sbinput.set_time(start_time = 0, end_time = 1e8, log_scale = True, nsteps = 2000)
# The next line will change the settings from the previous line.
sbinput.set_time(start_time = 0, end_time = 5e7, log_scale = False, dt = 1e5)

# Set mass grid.
#   mass_grid: Possible values are 'Full_Isochrone', 'Isochrone_Large', 'Large', 'Small'

sbinput.set_mass_grid(mass_grid = 'Full_Isochrone')

# dt: Time step for printing out spectra in yrs.

sbinput.dt_print_spectra(dt = 1e6)

# Selects model atmosphere for low resolution spectra output.
#   atmo: Possible values 'Planck', 'Lejeune', 'Lejuene-Schmutz', 'Lejeune-Hiller', 'Pauldrach-Hiller'

sbinput.set_atmo(atmo = 'Pauldrach-Hiller')

# Sets metallicity for the optical high res spectra library.
#   Options are: Z = 0.001, 0.008, 0.02, 0.04

sbinput.high_res_Z(Z = 0.02)

# UV_line_lib: [str] Choice of the UV spectral library, 'Solar' or 'LMC/SMC'

sbinput.uv_line_lib(lib = 'Solar')

# RGS microturbulance and abundance.
#   Possible values:
#     0 <= RGS_turb_vel <= 6
#     RGS_abund = True/False

sbinput.rgs_turb_abund(RGS_turb_vel = 3, RGS_abund = True)

# writeinput: Controls whether input parameters are written to a
#             textfile, separate from the model_name.input file 
#             written for Starburst99.
#                                
# silent:     Surpresses all file output. Data only kept in out_data 
#             object that can be used in python. Use this if you
#             only want to work with the output data in an active 
#             python session.
#             Overrides all other output choices. 
#                    
# Output files controled by True/False. By default all are True, except HRD (3).
#
# OUTPUT FILES:
#
#    1 model_name.QUANTA
#    2 model_name.SNR
#    3 model_name.HRD
#    4 model_name.POWER
#    5 model_name.SP
#    6 model_name.YIELDS
#    7 model_name.SPECTRUM
#    8 model_name.LINE
#    9 model_name.COLOR
#   10 model_name.WIDTH
#   11 model_name.FEATURES
#   12 model_name.OVI
#   13 model_name.HIRES
#   14 model_name.WRLINES
#   15 model_name.IFASPEC

sbinput.set_output(writeinput = False, 
                   silent = False, 
                   QUANTA = True, 
                   SNR = True, 
                   HRD = False, 
                   POWER = True, 
                   SP = True, 
                   YIELDS = True, 
                   SPECTRUM = True, 
                   LINE = True, 
                   COLOR = True, 
                   WIDTH = True, 
                   FEATURES = True, 
                   OVI = True, 
                   HIRES = True, 
                   WRLINES = True, 
                   IFASPEC = True)


sbinput.run_starburst()

