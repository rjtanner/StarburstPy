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

datapaths = sb.indata
datapaths.output_dir = './sboutput'
datapaths.SB99_dir = '../../starburst99'

model_name = 'happy'

sbinput = sb.SbInput(model_name, rank = 3)

sbinput.set_star_formation(fixed_mass = False, SFR = 1.0)
sbinput.set_star_formation(fixed_mass = True, total_mass = 1e6)

sbinput.set_IMF(IMF = 'Salpeter')
sbinput.set_IMF(IMF = 'Kroupa')
sbinput.set_IMF(IMF_exp = [1.3,2.3], IMF_bound = [0.1, 0.5, 120.0], mass_lo = 0.1, mass_hi = 120.0)

sbinput.set_SN_cut_off(lo_mass = 8.0, hi_mass = 120.0)

sbinput.set_evolve_track(track = 'Genevav40', Z = 0.002)

sbinput.set_wind_model(wind = 'Evolution')

sbinput.set_time(start_time = 0, end_time = 1e8, log_scale = True, nsteps = 2000)
sbinput.set_time(start_time = 0, end_time = 5e7, log_scale = False, dt = 1e5)

sbinput.set_mass_grid(mass_grid = 'Full_Isochrone')

sbinput.dt_print_spectra(dt = 1e6)

sbinput.set_atmo(atmo = 'Pauldrach-Hiller')

sbinput.high_res_Z(Z = 0.02)

sbinput.uv_line_lib(lib = 'Solar')

sbinput.rgs_turb_abund(RGS_turb_vel = 3, RGS_abund = True)

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

