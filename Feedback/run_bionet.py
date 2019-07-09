"""Simulates an example network of 450 cell receiving two kinds of exernal input as defined in the configuration file"""
import sys
from bmtk.simulator import bionet
from bmtk.simulator.bionet.modules.sim_module import SimulatorMod
from bmtk.utils.reports.spike_trains import SpikeTrains
from bmtk.utils.reports.spike_trains import PoissonSpikeGenerator

import numpy as np
import matplotlib.pyplot as plt
from neuron import h


firing_rate_threshold = 5.5


class FeedbackLoop(SimulatorMod):
    def __init__(self):
        self._synapses = {}
        self._netcon = None
        self._spike_events = None

        self._high_level_neurons = []
        self._low_level_neurons = []

        self._synapses = {}
        self._netcons = {}

        self._spike_records = {}

    def _set_spike_detector(self, sim):
        for gid in self._low_level_neurons:
            tvec = h.Vector()
            nc = self._netcons[gid]
            nc.record(tvec)
            self._spike_records[gid] = tvec

    def _activate_hln(self, sim, block_interval, firing_rate):
        next_block_tstart = (block_interval[1] + 1) * sim.dt  # The next time-step
        next_block_tstop = next_block_tstart + sim.nsteps_block*sim.dt  # The time-step when the next block ends
        # TODO: See if we've reached the end of the simulation

        # This is where you can use the firing-rate of the low-level neurons to generate a set of spike times for the
        # next block. For this case I'm just setting the high-level neuron to start bursting
        if firing_rate > firing_rate_threshold:
            self._spike_events = np.linspace(next_block_tstart, next_block_tstop, 500)
            #psg = PoissonSpikeGenerator()
            #psg.add(node_ids=[0], firing_rate=firing_rate, times=(next_block_tstart, next_block_tstop))
            #self._spike_events = psg.get_times(0)

            for gid in self._high_level_neurons:
                nc = self._netcons[gid]
                for t in self._spike_events:
                    nc.event(t)

    def initialize(self, sim):
        # Attach a NetCon/synapse on the high-level neuron(s) soma. We can use the NetCon.event(time) method to send
        # a spike to the synapse. Which, is a spike occurs, the high-level neuron will inhibit the low-level neuron.
        spikes = h.Vector([])  # start off with empty input
        vec_stim = h.VecStim()
        vec_stim.play(spikes)
        self._high_level_neurons = list(sim.net.get_node_set('high_level_neurons').gids())
        for gid in self._high_level_neurons:
            cell = sim.net.get_cell_gid(gid)

            # Create synapse
            # These values will determine how the high-level neuron behaves with the input
            new_syn = h.Exp2Syn(0.5, cell.hobj.soma[0])
            new_syn.e = 0.0
            new_syn.tau1 = 1.0
            new_syn.tau2 = 3.0

            nc = h.NetCon(vec_stim, new_syn)
            nc.weight[0] = 0.001
            nc.delay = 1.0

            self._synapses[gid] = new_syn
            self._netcons[gid] = nc

        # Attach another netcon to the low-level neuron(s) that will record
        self._low_level_neurons = list(sim.net.get_node_set('low_level_neurons').gids())
        for gid in self._low_level_neurons:
            cell = sim.net.get_cell_gid(gid)
            # NOTE: If you set the segment to 0.5 it will interfere with the record_spikes module. Neuron won't let
            # us record from the same place multiple times, so for now just set to 0.0. TODO: read and save spikes in biosimulator.
            nc = h.NetCon(cell.hobj.soma[0](0.0)._ref_v, None, sec=cell.hobj.soma[0])
            self._netcons[gid] = nc

        self._set_spike_detector(sim)

    def step(self, sim, tstep):
        pass

    def block(self, sim, block_interval):
        """This function is called every n steps during the simulation, as set in the config.json file (run/nsteps_block).

        We can use this to get the firing rate of neuron 1 during the last block and use it to calculate
        """

        # Calculate the avg number of spikes per neuron
        block_length = sim.nsteps_block*sim.dt/1000.0  #  time lenght of previous block of simulation TODO: precalcualte
        n_gids = 0
        n_spikes = 0
        for gid, tvec in self._spike_records.items():
            n_gids += 1
            n_spikes += len(list(tvec))  # use tvec generator. Calling this deletes the values in tvec

        avg_spikes = n_spikes/float(n_gids)

        # calculate the firing rate the the low-level neuron(s)
        fr = avg_spikes/float(block_length)
        #if fr > firing_rate_threshold:
        #    self._activate_hln(sim, block_interval, fr)

        # set the activity of high-level neuron
        self._activate_hln(sim, block_interval, fr)

        # NEURON requires resetting NetCon.record() every time we read the tvec.
        self._set_spike_detector(sim)

    def finalize(self, sim):
        pass


def run(config_file):
    conf = bionet.Config.from_json(config_file, validate=True)
    conf.build_env()

    fbmod = FeedbackLoop()

    graph = bionet.BioNetwork.from_config(conf)
    sim = bionet.BioSimulator.from_config(conf, network=graph)
    sim.add_mod(fbmod)  # Attach the above module to the simulator.
    sim.run()

    # Plotting
    spike_trains = SpikeTrains.from_sonata('output/spikes.h5')
    for gid in [0, 1]:
        spikes = np.zeros(sim.n_steps, dtype=np.int)
        if len(spike_trains.get_times(gid)) > 0:
            spikes[(spike_trains.get_times(gid)/sim.dt).astype(np.int)] = 1
        window_size = 10000
        window = np.zeros(window_size)
        window[0:window_size] = 1
        plt.plot(np.convolve(spikes, window), label=gid)
    plt.legend()
    plt.show()

    bionet.nrn.quit_execution()




if __name__ == '__main__':
    if __file__ != sys.argv[-1]:
        run(sys.argv[-1])
    else:
        run('config.json')
