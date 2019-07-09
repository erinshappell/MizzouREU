from bmtk.analyzer.cell_vars import plot_report
from scipy import signal

#plot_report(config_file='simulation_config.json')

import h5py
import matplotlib.pyplot as plt
import numpy as np

t = np.arange(0,1000,0.1)

### Plot Current Data ###

f = h5py.File('output/syn_report.h5')
i_data = [l[0] for l in list(f['data'])]
i_data = np.abs(i_data)

plt.figure()
plt.plot(t,i_data)
plt.xlabel('Time [ms]')
plt.ylabel('Synaptic Current [nA]')

### Plot Voltage Data ###

f2 = h5py.File('output/cell_vars.h5')
v_data = [l[0] for l in list(f2['data'])]

plt.figure()
plt.plot(t,v_data)
plt.xlabel('Time [ms]')
plt.ylabel('Membrane Voltage [mV]')

### Plot Voltage Spikes ###

#v_spikes = signal.find_peaks(v_data)
#plt.figure()
#plt.plot(v_spikes[0],np.ones(len(v_spikes[0])),'x')
#plt.xlim([0,10000])
#plt.xlabel('Sample #')
#plt.ylabel('Spikes')

### Plot both in same figure--show corresponding current and voltage spikes ###
fig, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('Time (t) [ms]')
ax1.set_ylabel('Synaptic Current (I) [nA]', color=color)
ax1.plot(t, i_data, color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
color = 'tab:blue'
ax2.set_ylabel('Membrane Voltage (V) [mV]', color=color)  # we already handled the x-label with ax1
ax2.plot(t, v_data, color=color)
ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()





 
