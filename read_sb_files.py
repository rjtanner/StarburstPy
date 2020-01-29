# -*- coding: utf-8 -*-
"""
    A script to read in data from a model run previously. Run a model using 
    'run_sbpy.py' or your own. Makes three plots.
    
    Author: Ryan Tanner
    email: ryan.tanner@nasa.gov
"""

import StarburstPy as sb
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['font.size'] = 18

datapaths = sb.indata
datapaths.output_dir = '../sboutput'

model_name = 'happy'
output_data = sb.out_data(model_name)
output_data.data, output_data.headers = sb.read_output_data(model_name = model_name, output_dir = datapaths.output_dir)

#fig = output_data.plot_power()

data = output_data.data

x = data['power']['Time']/1e6
y1 = 10**data['power']['Power_All']
y2 = 10**data['power']['Power_OB']
y3 = 10**data['power']['Power_RSG']
y4 = 10**data['power']['Power_LBV']
y5 = 10**data['power']['Power_WR']

fig = plt.figure(figsize = (12,8))
ax1 = fig.add_subplot(111)
ax1.semilogy(x, y1, label = 'Total Power',linewidth=10)
ax1.semilogy(x, y2, label = 'Power from OB Stars',linewidth=5)
ax1.semilogy(x, y3, label = 'Power from RSG Stars',linewidth=5)
ax1.semilogy(x, y4, label = 'Power from LBV Stars',linewidth=5)
ax1.semilogy(x, y5, label = 'Power from WR Stars',linewidth=5)
plt.xlabel('Time (Myrs)') 
plt.ylabel('All Power ergs/s')
plt.ylim(1e30,1e41)
plt.legend(loc='lower left')

plt.show()
input("Press Enter to continue...")


x = data['snr']['Time']/1e6
y2 = 10**data['snr']['Type_IB_Power']
y3 = 10**data['snr']['Stars_SN_Power']

fig = plt.figure(figsize = (12,8))
ax1 = fig.add_subplot(111)
ax1.plot(x, y2, label = 'Power from Type Ib SN',linewidth=5)
ax1.plot(x, y3, label = 'Power from Stellar SN',linewidth=5)

plt.xlabel('Time (Myrs)') 
plt.ylabel('All Power ergs/s')
plt.ylim(1e40,4e40)
plt.legend(loc='upper right')

plt.show()
input("Press Enter to continue...")


x = data['chemical_yields']['Time']/1e6
y1 = 10**data['chemical_yields']['H']
y2 = 10**data['chemical_yields']['He']
y3 = 10**data['chemical_yields']['C']
y4 = 10**data['chemical_yields']['N']
y5 = 10**data['chemical_yields']['O']
y6 = 10**data['chemical_yields']['Mg']
y7 = 10**data['chemical_yields']['Si']
y8 = 10**data['chemical_yields']['S']
y9 = 10**data['chemical_yields']['Fe']

fig = plt.figure(figsize = (12,8))
ax1 = fig.add_subplot(111)
ax1.semilogy(x, y1, label = 'H',linewidth=5)
ax1.semilogy(x, y2, label = 'He',linewidth=5)
ax1.semilogy(x, y3, label = 'C',linewidth=5)
ax1.semilogy(x, y4, label = 'N',linewidth=5)
ax1.semilogy(x, y5, label = 'O',linewidth=5)
ax1.semilogy(x, y6, label = 'Mg',linewidth=5)
ax1.semilogy(x, y7, label = 'Si',linewidth=5)
ax1.semilogy(x, y8, label = 'S',linewidth=5)
ax1.semilogy(x, y9, label = 'Fe',linewidth=5)
plt.xlabel('Time (Myrs)') 
plt.ylabel('Mass Loss (Msun/yr)')
#plt.ylim(1e30,1e41)
plt.legend(loc='lower left')

plt.show()
input("Press Enter to continue...")