# -*- coding: utf-8 -*-
"""
    A script to read in the data from models run in 'run_multi_models.py'.
    Picks off the total power of the starburst and plots the power vs. time.
    
    ****Need to run 'run_multi_models.py' FIRST****
    
    Author: Ryan Tanner
    email: ryan.tanner@nasa.gov

"""


import StarburstPy as sb
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['font.size'] = 18

datapaths = sb.indata
datapaths.output_dir = '../sboutput/sfr'
model_name_base = 'happy'

sfrrange = np.logspace(1,2,num=10)
y = []

for sfr in sfrrange:
    model_name = model_name_base+'{0}'.format(int(sfr))
    output_data = sb.out_data(model_name)
    output_data.data, output_data.headers = sb.read_output_data(model_name = model_name, output_dir = datapaths.output_dir)
    data = output_data.data
    x = data['power']['Time']/1e6
    y.append(10**data['power']['Power_All'])
    
fig = plt.figure(figsize = (12,8))
ax1 = fig.add_subplot(111)

for i in range(10):
    ax1.semilogy(x, y[i], label = 'SFR {0}'.format(int(sfrrange[i])),linewidth=5)

plt.xlabel('Time (Myrs)') 
plt.ylabel('Total Power ergs/s')
plt.ylim(1e40,1e43)
plt.legend(loc='lower right')

plt.show()