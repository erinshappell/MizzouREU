import os
import numpy as np

###################################################################################################
# eus_rate() will generate a spike train corresponding to the filling/voiding of the bladder
# INPUTS: filename  -- name of the file to write bladder spikes to
#		  input_dir -- directory containing the file to be written to
#		  cells     -- array of cell indices
#		  hz		-- firing rates (high,low,high)
#		  start     -- start times for simulation
#		  end  		-- end times of simulation

def eus_rate(filename,input_dir,cells,hz,start,end):
    stims = {}   
    for i in np.arange(0,len(hz)):
        interval = 1000.0/hz[i]
        for cell in cells:
            spk = start[i]
            spike_ints = np.random.poisson(interval,100)
            for spike_int in spike_ints:
                if not stims.get(cell):
                    stims[cell] = []
                spk += spike_int
                if spk <= end[i]:
                    stims[cell].append(spk)
    
    with open(os.path.join(input_dir,filename),'w') as file:
        file.write('gid spike-times\n')
        for key, value in stims.items():
            file.write(str(key)+' '+','.join(str(x) for x in value)+'\n')   
    return

###################################################################################################
# pag_rate() will generate a spike train corresponding to the filling/voiding of the bladder
# INPUTS: filename  -- name of the file to write bladder spikes to
#		  input_dir -- directory containing the file to be written to
#		  cells     -- array of cell indices
#		  hz		-- firing rates (low,high,low)
#		  start     -- start times for simulation
#		  end  		-- end times of simulation

def pag_rate(filename,input_dir,cells,hz,start,end):
    stims = {}   
    for i in np.arange(0,len(hz)):
        interval = 1000.0/hz[i]
        for cell in cells:
            spk = start[i]
            spike_ints = np.random.poisson(interval,100)
            for spike_int in spike_ints:
                if not stims.get(cell):
                    stims[cell] = []
                spk += spike_int
                if spk <= end[i]:
                    stims[cell].append(spk)
    
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
	begin_void = 0
	end_void = 0
	v_t = v_0 + fill*spk*150 # Current volume of bladder
	# 150 is a scaling value needed for this simulation.

	# Collect volume data
	time = []
	volume = []

	for cell in cells:
		spk = start
		v_t = v_0 + fill*spk
		
		# Minimum number of spikes per interval size
		num_ints = int(duration/1000)
		
		# Filling
		while spk <= mid and v_t < max_v:
			interval = 10.0/v_t

			spike_ints = np.random.poisson(interval, num_ints)

			for spike_int in spike_ints:
				if not stims.get(cell):
					stims[cell] = []

				spk += spike_int

				# Update current volume of bladder
				vol = v_0 + fill*spk*150
				v_t = v_0 + fill*spk*150*v_t
				# 150 is a scaling value needed for this simulation.

				if spk <= duration:
					stims[cell].append(spk)

					# Save time and volume data for cell 0
					if cell == 0:
						time.append(spk)
						volume.append(vol)		

		# Update value for the running total of times at which filling ended
		# (and voiding began)
		begin_void += spk

		# Reset void flag (this will give us the running total of times
		# at which voiding ended)
		end_void_flag = 0

		# Voiding
		while spk <= duration:
			
			interval = 10.0/v_t
			spike_ints = np.random.poisson(interval, num_ints)
			
			for spike_int in spike_ints:
				if not stims.get(cell):
					stims[cell] = []

				spk += spike_int

				# Update current volume of bladder
				if v_t > v_0:
					v_t -= void*spk
				if v_t < v_0:
					v_t = v_0
					if end_void_flag == 0:
						end_void_flag = 1
						end_void += spk

				vol = v_t

				if spk <= duration:
					stims[cell].append(spk)

					# Save time and volume data for cell 0
					if cell == 0:
						time.append(spk)
						volume.append(vol)
	
	# Write spike data to file
	with open(os.path.join(input_dir,filename),'w') as file:
		file.write('gid spike-times\n')
		for key, value in stims.items():
			file.write(str(key)+' '+','.join(str(x) for x in value)+'\n')   
	
	begin_void /= len(cells) # Begin void time is average of all times for the 10 cells
	end_void   /= len(cells) # End void time is average of all times for the 10 cells

	return (begin_void,end_void,time,volume)

###################################################################################################

if __name__ == '__main__':
    
	# Change these for your needs	
	start    = 0 	 #ms
	mid 	 = 6000	 #ms
	duration = 10000 #ms
	cells = np.arange(0,250)

	# Create bladder afferent input spikes ----------------

	output_file = 'Blad_spikes.csv'
	input_dir = './'

	fill = 0.05 	 	# ml/min (Asselt et al. 2017)
	fill /= (1000 * 60) # Scale from ml/min to ml/ms
	void = 4.6 	 		# ml/min (Streng et al. 2002)
	void /= (1000 * 60) # Scale from ml/min to ml/ms
	max_v = 0.76 		# ml (Grill et al. 2019)

	# After bladder is full, how long should we wait before conscious activation of voiding?
	delay = 1000.0 # ms
	
	(begin_void,end_void,v_time,v_vol) = bladder_rate(output_file, input_dir, fill, void, max_v, cells, start, mid + delay, duration)

	# Create EUS afferent input spikes --------------------
	output_file = 'EUS_spikes.csv'
	input_dir = './'

	hz = [15.0,1.0,15.0]  # Using 1.0 Hz as basal firing rate of pudendal (EUS) afferent
					 	  # (Habler et al. 1993)
					 	  # Using high PAG firing rate of 15.0 Hz as high firing rate 
					 	  # (Blok et al. 2000)
	start = [0.0, begin_void - delay, end_void + delay]
	end = [begin_void - delay , end_void + delay, 10000.0]

	eus_rate(output_file, input_dir, cells, hz, start, end)

	# Create PAG afferent input spikes --------------------
	output_file = 'PAG_spikes.csv'
	input_dir = './'

	hz = [1.0,15.0,1.0] # Using basal firing rate of pudendal (EUS) afferent for PAG afferent
				   		# 1.0 Hz (Habler et al. 1993)
				   		# Using 15.0 Hz as high firing rate of PAG afferent
				   		# (Blok et al. 2000)
	start = [0.0, begin_void - delay, end_void + delay]
	end= [begin_void - delay , end_void + delay, 10000.0]
	pag_rate(output_file, input_dir, cells, hz, start, end)