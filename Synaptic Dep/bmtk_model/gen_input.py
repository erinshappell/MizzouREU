from bmtk.utils.spike_trains import SpikesGenerator
import numpy as np
import os

def abrupt_changing_rates(filename,input_dir,hz,start,end):
    #print('len(hz)=',len(hz))
    stims = {}  
    cell = 0
    for i in np.arange(0,len(hz)):
        interval = 1000.0/hz[i]
        spk = start[i]
        spike_ints = np.random.poisson(interval,100)
        for spike_int in spike_ints:
            if not stims.get(cell):
                stims[cell] = []
            spk = spk + spike_int
            if spk <= end[i]:
                stims[cell].append(spk)
    
    with open(os.path.join(input_dir,filename),'w') as file:
        file.write('gid spike-times\n')
        for key, value in stims.items():
            file.write(str(key)+' '+','.join(str(x) for x in value)+'\n')   
    return

def write_input_file(rate=50):
	sg = SpikesGenerator(nodes='network/input_nodes.h5', t_max=1.0)
	sg.set_rate(rate)
	sg.save_csv('./input/input_spikes.csv', in_ms=True)
	#output_file = 'input/input_spikes.csv'
	#input_dir = './'
	#hz = [15.0,30.0]
	#start = [0.0,500.0]
	#end = [500.0,1000.0]
	#abrupt_changing_rates(output_file, input_dir, hz, start, end)

if __name__ == '__main__':
    write_input_file()
