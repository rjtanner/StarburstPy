# -*- coding: utf-8 -*-
"""
    A script to run multiple Starburst99 models. Runs 10 models with continuous
    star formation with SFR from 1.0 to 100.0 spaced logarithmically.
    
    Author: Ryan Tanner
    email: ryan.tanner@nasa.gov
"""
import StarburstPy as sb
import numpy as np

datapaths = sb.indata
datapaths.output_dir = '../sboutput/sfr'
datapaths.SB99_dir = '../../starburst99'
model_name = 'happy'

sbinput = sb.SbInput('happy', chatter = 2)

sbinput.set_time(start_time = 0, end_time = 1e7, log_scale = False, dt = 1e5)
sbinput.set_output(out_type = 'original')

sfrrange = np.logspace(1,2,num=10)

for sfr in sfrrange:
    sbinput.set_model_name(model_name = model_name+'{0}'.format(int(sfr)))
    sbinput.set_star_formation(fixed_mass = False, SFR = sfr)
    print(model_name+'{0}'.format(int(sfr)))
    output_data = sbinput.run_starburst()