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
		
		# Filling (0.05 ml/min) (Asselt et al. 2017)
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

###################################################################################################
# bladder_rate() will generate a spike train corresponding to the filling/voiding of the bladder
# INPUTS: filename  -- name of the file to write bladder spikes to
#		  input_dir -- directory containing the file to be written to
#		  fill      -- fill rate for the bladder
#		  void      -- void rate for the bladder
#		  max_v     -- maximum volume of the bladder
#		  cells     -- array of cell indices
#		  start     -- start time for simulation
#		  mid       -- point where bladder transitions from filling to voiding
#		  duration  -- total time of simulation

def bladder_rate(filename,input_dir,fill,void,max_v,cells,start,mid,duration):
	
	stims = {}
	v_0 = 0.05 # Initial volume of bladder (needs checking)
	spk = start
	actual_mid = 0
	v_t = v_0 + fill*spk # Current volume of bladder

	for cell in cells:
		spk = start
		v_t = v_0 + fill*spk
		
		# Minimum number of spikes per interval size
		num_ints = int(duration/1000)
		
		# Filling
		while spk <= mid and v_t < max_v:
			interval = 10.0/v_t
			# 10.0 is a scaling value needed for this simulation
			spike_ints = np.random.poisson(interval, num_ints)
			if cell == 0:
				print('interval is %d' %interval)
				print(spike_ints)

			for spike_int in spike_ints:
				if not stims.get(cell):
					stims[cell] = []
				spk += spike_int
				if spk <= duration:
					stims[cell].append(spk)

			# Update current volume of bladder
			v_t = v_0 + fill*spk*v_t*150 
			# 150 is a scaling value needed for this simulation.

			if cell == 0:
				print('spk is %d' %spk)
				print('volume is %.2f' %v_t)

		v_end = v_t
		actual_mid += spk

		if cell == 0:
			print('end of filling...')
			print('v_end is %.2f' %v_end)

		# Voiding
		while spk <= duration:
			interval = 10.0/v_t
			spike_ints = np.random.poisson(interval, num_ints)

			if cell == 0:
				print('interval is %d' %interval)
				print(spike_ints)
			
			for spike_int in spike_ints:

				if not stims.get(cell):
					stims[cell] = []
				spk += spike_int
				if spk <= duration:
					stims[cell].append(spk)

			# Update current volume of bladder
			if v_t > v_0:
				v_t -= void*spk
			if v_t < v_0:
				v_t = v_0
		
			if cell == 0:
				print('spk is %d' %spk)
				print('volume is %.2f' %v_t)
			
	
	# Write spike data to file
	with open(os.path.join(input_dir,filename),'w') as file:
		file.write('gid spike-times\n')
		for key, value in stims.items():
			file.write(str(key)+' '+','.join(str(x) for x in value)+'\n')   
	
	actual_mid /= len(cells) # Actual mid is average of all mids for the 10 cells
	return actual_mid

###################################################################################################

if __name__ == '__main__':
    
	# Change these for your needs	
	start    = 0 	 #ms
	mid 	 = 6000	 #ms
	duration = 10000 #ms
	cells = [0,1,2,3,4,5,6,7,8,9]

	# Create bladder afferent input spikes ----------

	output_file = 'Blad_spikes.csv'
	input_dir = './'

	fill = 0.05 	 	# ml/min (Asselt et al. 2017)
	fill /= (1000 * 60) # Scale from ml/min to ml/ms
	void = 4.6 	 		# ml/min (Streng et al. 2002)
	void /= (1000 * 60) # Scale from ml/min to ml/ms
	max_v = 0.76 		# ml (Grill et al. 2019)

	# After bladder is full, how long should we wait before conscious activation of voiding?
	delay = 500.0 # ms
	
	actual_mid = bladder_rate(output_file, input_dir, fill, void, max_v, cells, start, mid + delay, duration)

	# Create EUS afferent input spikes --------------
	output_file = 'EUS_spikes.csv'
	input_dir = './'

	hz = [15.0,1.0]  # Using 1.0 Hz as basal firing rate of pudendal (EUS) afferent
					 # (Habler et al. 1993)
					 # Using high PAG firing rate of 15.0 Hz as high firing rate 
					 # (Blok et al. 2000)
	start = [0.0, actual_mid]
	end= [actual_mid, 10000.0]

	abrupt_changing_rates(output_file, input_dir, hz, start, end, cells)

	# Create PAG afferent input spikes --------------
	output_file = 'PAG_spikes.csv'
	input_dir = './'

	hz = [1.0,15.0] # Using basal firing rate of pudendal (EUS) afferent for PAG afferent
				   # 1.0 Hz (Habler et al. 1993)
				   # Using 15.0 Hz as high firng rate of PAG afferent
				   # (Blok et al. 2000)
	start = [0.0, actual_mid - delay]
	end= [actual_mid - delay, 10000.0]
	abrupt_changing_rates(output_file, input_dir, hz, start, end, cells)

	
