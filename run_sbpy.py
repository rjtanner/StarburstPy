# -*- coding: utf-8 -*-

import StarburstPy as sb

"""
    Author: Ryan Tanner
    email: ryan.tanner@nasa.gov
"""

# Sets path links to current directory, output directory, and the directory where
# the Starburst99 fortran code is compiled. Defaults set to current directory.

datapaths = sb.indata
datapaths.output_dir = './sboutput'
datapaths.SB99_dir = '../../starburst99'

# A model name must be passed in to make the input object.

model_name = 'happy'

# The input object "sbinput" contains all the input parameters and methods for 
# setting all parameters, and a method for running Starburst99.
# If just the model name is given then the default parameters are set to the
# defaults from the Starburst99 webpage.

sbinput = sb.SbInput(model_name)

# Method for running Starburst99 based on input parameters that have been set.

output_data = sbinput.run_starburst()

# Data from Starburst99 is in
# output_data.data
#
# Headers from the output files are found in
# output_data.headers
#
# See run_sbpy_all_options.py for examples of all possible input parameters.
