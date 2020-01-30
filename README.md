# Starburstpy

StarburstPy is a python wrapper for Starburst99 (Leitherer et al. (1999), Leitherer et al. (2014)) in python.
    
The Starburst99 Fortan code must be downloaded and compiled before this will work. See http://www.stsci.edu/science/starburst99/docs/default.htm for the Starburst99 code and instructions on how to compile it.
    
I do not maintain the Starburst99 code so I cannot answer (most) technical questions about it. Questions about Starburst99 should be directed to Claus Leitherer.

Example scripts can be found at https://github.com/rjtanner/StarburstPy. Code can also be found at https://pypi.org/project/StarburstPy/.

StarburstPy can also be used to read in Starburst99 data files generated previously. See 'read_sb_files.py' for an example.

Download from github or install using 

pip install StarburstPy

Questions about StarburstPy should be directed to me, the author. 
    
    Author: Ryan Tanner
    email: ryan.tanner@nasa.gov
    url: https://github.com/rjtanner/StarburstPy
    
    
v0.9.1

Output data is stored in its own class. It contains two dictionaries. One with the data and one with file header information.

Output files 'irfeature', 'sptyp2', and 'wrlines' are not supported yet. Will include in v0.9.2.
