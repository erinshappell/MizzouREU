import os
import numpy as np

from bmtk.builder.networks import NetworkBuilder


net = NetworkBuilder("internal")

# Create a inhibitory cell (the high-level neuron) gid = 0
net.add_nodes(N=1, level='high', ei='i',
              model_type='biophysical',
              model_template='ctdb:Biophys1.hoc',
              model_processing='aibs_perisomatic',
              dynamics_params='472912177_fit.json',
              morphology='Pvalb-IRES-Cre_Ai14_IVSCC_-176847.04.02.01_470522102_m.swc')

# Create an excitatory cell (the low-level neuron) gid = 1.
net.add_nodes(N=1, level='low', ei='e',
              model_type='biophysical',
              model_template='ctdb:Biophys1.hoc',
              model_processing='aibs_perisomatic',
              dynamics_params='473863035_fit.json',
              morphology='Nr5a1-Cre_Ai14_IVSCC_-169250.03.02.01_471087815_m.swc')

# Create synapse from the inhibitory cell to the excitatory cell (0 --> 1)
net.add_edges(source={'level': 'high'},
              target={'level': 'low'},
              connection_rule=5,
              syn_weight=0.8,
              target_sections=['somatic'],
              distance_range=[0.0, 1e+20],
              delay=1.0,
              dynamics_params='GABA_InhToExc.json',
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
              syn_weight=10.0,
              # weight_function='wmax',
              distance_range=[0.0, 10.0],
              target_sections=['somatic'],
              delay=1.0,
              dynamics_params='AMPA_ExcToExc.json',
              model_template='exp2syn')

inputs.build()
inputs.save(output_dir='network')
