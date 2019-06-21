from bmtk.analyzer.cell_vars import plot_report
from bmtk.analyzer.spike_trains import raster_plot
import pandas as pd
import matplotlib.pyplot as plt
from bmtk.utils.cell_vars import CellVarsFile
import numpy as np
from bmtk.analyzer.cell_vars import _get_cell_report
#plot_report(config_file="simulation_config.json")
import pdb


numBladaff = 10
numEUSaff = 10
numPAGaff = 10
numIND = 10
numHypo = 10
numINmplus = 10
numINmminus = 10
numPGN = 10
numFB = 10
numIMG = 10
numMPG = 10
numEUSmn = 10
numBladmn = 10

config_file = "simulation_config.json"
report_name = None
report_name, report_file = _get_cell_report(config_file, report_name)


var_report = CellVarsFile(report_file)
time_steps = var_report.time_trace





# Plot spike raster

Blad_gids = np.arange(0,numBladaff)
EUS_gids = Blad_gids + numBladaff
PAG_gids = EUS_gids + numEUSaff
IND_gids = PAG_gids + numPAGaff
Hypo_gids = IND_gids + numIND
INmplus_gids = Hypo_gids + numHypo
INmminus_gids = INmplus_gids + numINmplus
PGN_gids = INmminus_gids + numINmminus
FB_gids = PGN_gids + numPGN
IMG_gids = FB_gids + numFB
MPG_gids = IMG_gids + numIMG
EUSmn_gids = MPG_gids + numMPG
Bladmn_gids = EUSmn_gids + numEUSmn

df = pd.read_csv("output/spikes.csv",delimiter=' ')
rast = df.values

#plt.plot(rast[:,0],rast[:,1],'k.')

# Plot Bladder afferent, EUS afferent, PAG afferent, IND, and Hypo on one figure
plt.figure()
Bladspkt = rast[np.in1d(rast[:,1],Blad_gids),:]
plt.plot(Bladspkt[:,0],Bladspkt[:,1],'b.',label='Blad. aff')

EUSspkt = rast[np.in1d(rast[:,1],EUS_gids),:]
plt.plot(EUSspkt[:,0],EUSspkt[:,1],'r.',label='EUS aff')

PAGspkt = rast[np.in1d(rast[:,1],PAG_gids),:]
plt.plot(PAGspkt[:,0],PAGspkt[:,1],'m.',label='PAG aff')



plt.xlabel('Time (t) [ms]')
plt.legend()

# Plot INM+, INM-, IND, FB
plt.figure()
INMpspkt = rast[np.in1d(rast[:,1],INmplus_gids),:]
plt.plot(INMpspkt[:,0],INMpspkt[:,1]-40,'b.',label='INm+')

INMmspkt = rast[np.in1d(rast[:,1],INmminus_gids),:]
plt.plot(INMmspkt[:,0],INMmspkt[:,1]-40,'r.',label='INM-')

INDspkt = rast[np.in1d(rast[:,1],IND_gids),:]
plt.plot(INDspkt[:,0],INDspkt[:,1]-30,'g.',label='IND')

FBspkt = rast[np.in1d(rast[:,1],FB_gids),:]
plt.plot(FBspkt[:,0],FBspkt[:,1]-50,'k.',label='FB')

plt.xlabel('Time (t) [ms]')
plt.legend()

# Plot SPN(PGN), Hypo, MPG, IMG
plt.figure()
PGNspkt = rast[np.in1d(rast[:,1],PGN_gids),:]
plt.plot(PGNspkt[:,0],PGNspkt[:,1]-50,'b.',label='PGN')

IMGspkt = rast[np.in1d(rast[:,1],IMG_gids),:]
plt.plot(IMGspkt[:,0],IMGspkt[:,1]-80,'r.',label='IMG')

MPGspkt = rast[np.in1d(rast[:,1],MPG_gids),:]
plt.plot(MPGspkt[:,0],MPGspkt[:,1]-70,'c.',label='MPG')

Hypospkt = rast[np.in1d(rast[:,1],Hypo_gids),:]
plt.plot(Hypospkt[:,0],Hypospkt[:,1]-40,'m.',label='Hypo')

plt.xlabel('Time (t) [ms]')
plt.legend()

# Motor neurons
plt.figure()
EUSmnspkt = rast[np.in1d(rast[:,1],EUSmn_gids),:]
plt.plot(EUSmnspkt[:,0],EUSmnspkt[:,1],'r.',label='EUSmn')

Bladmnspkt = rast[np.in1d(rast[:,1],Bladmn_gids),:]
plt.plot(Bladmnspkt[:,0],Bladmnspkt[:,1],'g.',label='Bladmn')

plt.xlabel('Time (t) [ms]')
plt.legend()




#Plot averaged frequency
#plt.figure()
#plt.hist(Bladspkt[:,0],10)

#plt.figure()
#plt.hist(Bladmnspkt[:,0],10)

plt.show()
#raster_plot(cells_file="network/inputPUD_hco_net_edges.h5", cell_models_file="network/inputPUD_hco_net_edge_types.csv", spikes_file="output/spikes.h5")
