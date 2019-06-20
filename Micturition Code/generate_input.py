import os
import numpy as np

def constant_rate_input(filename,input_dir,hz,start,duration,cells):
    
    stims = {}   
    interval = 1000.0/hz
    for cell in cells:
        spk = start
        spike_ints = np.random.poisson(interval,100)
        for spike_int in spike_ints:
            if not stims.get(cell):
                stims[cell] = []
            if spk <= duration:
                spk = spk + spike_int
                stims[cell].append(spk)
    
    with open(os.path.join(input_dir,filename),'w') as file:
        file.write('gid spike-times\n')
        for key, value in stims.items():
            file.write(str(key)+' '+','.join(str(x) for x in value)+'\n')   
    return


def abrupt_changing_rates(filename,input_dir,hz,start,end,cells):
    #print('len(hz)=',len(hz))
    stims = {}   
    for i in np.arange(0,len(hz)):
        interval = 1000.0/hz[i]
        for cell in cells:
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

def ramp_rate_input(filename,input_dir,hz,start,duration,cells,divs):
    
	stims = {}
	scale_step = (hz[1]-hz[0])/divs
	for cell in cells:
		spk = start
		
		# Minimum number of spikes per interval size
		num_ints = int(duration/1000)
		
		# Scaling factor for interval size
		int_scale = hz[0]

		while spk <= duration:
			interval = 1000.0/int_scale
			spike_ints = np.random.poisson(interval, num_ints + int(interval/1000))
			
			# Max interval scaling factor is the max frequency
			if int_scale > hz[1]:
				int_scale = hz[1]

			for spike_int in spike_ints:
				if not stims.get(cell):
					stims[cell] = []
				spk += spike_int
				if spk <= duration:
					stims[cell].append(spk)

			# Update interval scaling factor
			int_scale += scale_step
  
	# Write spike data to file
	with open(os.path.join(input_dir,filename),'w') as file:
		file.write('gid spike-times\n')
		for key, value in stims.items():
			file.write(str(key)+' '+','.join(str(x) for x in value)+'\n')   
		return

def ramp_rate_updown(filename,input_dir,hz,start,duration,cells,divs):
	
	stims = {}
	scale_step = (hz[1]-hz[0])/divs
	mid = int(duration/2)
	for cell in cells:
		spk = start
		
		# Minimum number of spikes per interval size
		num_ints = int(duration/1000)
		
		# Scaling factor for interval size
		int_scale = hz[0]

		
		# Filling (0.05 ml/min) (Asselt et al. 2006)
		while spk <= mid:
			interval = 1000.0/int_scale
			spike_ints = np.random.poisson(interval, num_ints + int(interval/1000))

			for spike_int in spike_ints:
				if not stims.get(cell):
					stims[cell] = []
				spk += spike_int
				if spk <= duration:
					stims[cell].append(spk)

			# Update interval scaling factor
			if int_scale < hz[1]:
				int_scale += scale_step

			# Max interval scaling factor is the max frequency
			if int_scale > hz[1]:
				int_scale = hz[1]


		int_scale -= scale_step
		num_ints = int(num_ints/10)

		# Voiding (5 ml/min) (Streng et al. 2002)
		while spk <= duration:
			interval = 1000.0/int_scale
			spike_ints = np.random.poisson(interval, num_ints)
			
			# Max interval scaling factor is the max frequency
			if int_scale > hz[1]:
				int_scale = hz[1]

			for spike_int in spike_ints:
				if not stims.get(cell):
					stims[cell] = []
				spk += spike_int
				if spk <= duration:
					stims[cell].append(spk)

			# Update interval scaling factor
			if int_scale > hz[0]:
				int_scale -= scale_step

			# Min interval scaling factor is the min frequency
			if int_scale < hz[0]:
				int_scale = hz[0]
  
	# Write spike data to file
	with open(os.path.join(input_dir,filename),'w') as file:
		file.write('gid spike-times\n')
		for key, value in stims.items():
			file.write(str(key)+' '+','.join(str(x) for x in value)+'\n')   
		return

if __name__ == '__main__':
    
	#Change these for your needs

	# Create bladder afferent input spikes ----------
	start = 0 #ms
	mid = 5000 #ms
	duration = 10000 #ms

	output_file = 'Blad_spikes.csv'
	input_dir = './'

	hz = [4.0,16.0] # Starting with 4 Hz then moving up 16 Hz as bladder fills 
			# (Grill et al. 2016)
	cells = [0,1,2,3,4,5,6,7,8,9]
	# Number of divisions = number of different rates between/including
	# hz(min) and hz(max)
	divs = 4
	#ramp_rate_input(output_file, input_dir, hz, start, duration, cells, divs)
	ramp_rate_updown(output_file, input_dir, hz, start, duration, cells, divs)

	# Create EUS afferent input spikes --------------
	output_file = 'EUS_spikes.csv'
	input_dir = './'

	hz = [16.0,4.0]
	start = [0.0,5000.0]
	end= [5000.0,10000.0]
	abrupt_changing_rates(output_file, input_dir, hz, start, end, cells)

	# Create PAG afferent input spikes --------------
	output_file = 'PAG_spikes.csv'
	input_dir = './'

	hz = [4.0,16.0]
	start = [0.0,5000.0]
	end= [5000.0,10000.0]
	abrupt_changing_rates(output_file, input_dir, hz, start, end, cells)

	
