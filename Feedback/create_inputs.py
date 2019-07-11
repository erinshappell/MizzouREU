import os
from bmtk.utils.reports.spike_trains import PoissonSpikeGenerator
from bmtk.utils.reports.spike_trains import SpikeTrains
import matplotlib.pyplot as plt
import numpy as np

if os.path.exists('inputs/ramping_spikes.h5'):
    os.remove('inputs/ramping_spikes.h5')

# Create PGN firing rate
time_range = np.linspace(0.0, 10.0, 20000)
psg = PoissonSpikeGenerator(mode='w')
psg.add(node_ids=[0], firing_rate=np.linspace(1.0, 10.0, 30000), population='external', times=time_range)
psg.to_sonata('inputs/ramping_spikes.h5')

#st = SpikeTrains.from_sonata('inputs/ramping_spikes.h5')
#print(st.get_times(0))
#spikes = np.zeros(100001, dtype=np.int)
#spikes[(st.get_times(0) / 0.1).astype(np.int)] = 1
#plt.plot(spikes)
#plt.show()