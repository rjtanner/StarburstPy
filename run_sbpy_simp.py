# -*- coding: utf-8 -*-

import StarburstPy as sb

"""
    A simple script with the major options.

    Author: Ryan Tanner
    email: ryan.tanner@nasa.gov

"""

datapaths = sb.indata

datapaths.output_dir = './sboutput'
datapaths.SB99_dir = '../../starburst99'

model_name = 'happy'

sbinput = sb.SbInput(model_name)

#sbinput.set_star_formation(fixed_mass = False, SFR = 15.0)
sbinput.set_star_formation(fixed_mass = True, total_mass = 1e6)

sbinput.set_IMF(IMF = 'Kroupa')

sbinput.set_evolve_track(track = 'Genevav40', Z = 0.002)

#sbinput.set_time(start_time = 0, end_time = 1e7, log_scale = True, nsteps = 1000)
sbinput.set_time(start_time = 0, end_time = 1e7, log_scale = False, dt = 1e5)


# Note: The writeinput option will make a txt file *separate* from the standard
#       Starburst99 input file.

sbinput.set_output(writeinput = False, 
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
