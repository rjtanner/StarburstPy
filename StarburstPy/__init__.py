# -*- coding: utf-8 -*-
"""
    StarburstPy is a python wrapper for Starburst99 (Leitherer et al. (1999), 
    Leitherer et al. (2014)) in python.
    
    The Starburst99 Fortan code must be downloaded and compiled before this 
    will work. See http://www.stsci.edu/science/starburst99/docs/default.htm
    for the Starburst99 code and instructions on how to compile it.
    
    I do not maintain the Starburst99 code so I cannot answer (most) technical 
    questions about it. Questions about Starburst99 should be directed to 
    Claus Leitherer.

    Questions about StarburstPy should be directed to me, the author.    
    
    Author: Ryan Tanner
    email: ryan.tanner@nasa.gov

"""

__all__ = ['starburstpy', 'utils']

from .utils import *
from .utils.sb_stdout import sb_messages
from .utils.file_io import *
from .starburstpy import SbInput, out_data, run_starburst

sb_mess = sb_messages(rank = 2)

indata = _file_paths()