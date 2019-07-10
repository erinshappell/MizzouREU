import os
import numpy as np

from bmtk.builder.networks import NetworkBuilder


net = NetworkBuilder("internal")

# Represents the bladder afferent
# Create a excitatory cell (the high-level neuron) gid = 0
net.add_nodes(N=1, level='high', ei='i',
              model_type='biophysical',
              model_template='ctdb:Biophys1.hoc',
              model_processing='aibs_perisomatic',
              dynamics_params='472912177_fit.json',
              morphology='Pvalb-IRES-Cre_Ai14_IVSCC_-176847.04.02.01_470522102_m.swc')

# Represents the PGN
# Create an excitatory cell (the low-level neuron) gid = 1.
net.add_nodes(N=1, level='low', ei='e',
              model_type='biophysical',
              model_template='ctdb:Biophys1.hoc',
              model_processing='aibs_perisomatic',
              dynamics_params='473863035_fit.json',
              morphology='Nr5a1-Cre_Ai14_IVSCC_-169250.03.02.01_471087815_m.swc')

# Create synapse from the high-level cell to the low-level cell (0 --> 1)
net.add_edges(source={'level': 'high'},
              target={'level': 'low'},
              connection_rule=5,
              syn_weight=12.0e-03,
              target_sections=['somatic'],
              distance_range=[0.0, 300.0],
              delay=1.0,
              dynamics_params='AMPA_ExcToExc.json',
              model_template='exp2syn')

# build and save the network.
net.build()
net.save(output_dir='network')

# An external network with 1 (excitatory virtual) cell that will connect to the internal excitatory cell to drive
# it with the spikes inputs
inputs = NetworkBuilder("external")
inputs.add_nodes(N=1,
              model_type='virtual')
inputs.add_edges(source=inputs.nodes(), target=net.nodes(level='low'),
              connection_rule=1,
              syn_weight=12.0E-03,
              # weight_function='wmax',
              distance_range=[0.0, 300.0],
              target_sections=['somatic'],
              delay=1.0,
              dynamics_params='AMPA_ExcToExc.json',
              model_template='exp2syn')

inputs.build()
inputs.save(output_dir='network')
