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
		actual_mid += spk

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
	
	actual_mid /= len(cells) # Actual mid is average of all mids for the 10 cells
	return (actual_mid,time,volume)

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
	
	(actual_mid,v_time,v_vol) = bladder_rate(output_file, input_dir, fill, void, max_v, cells, start, mid + delay, duration)

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

	# Plot volume -----------------------------------
	import matplotlib.pyplot as plt 

	v_vol = np.array(v_vol)

	plt.figure()
	plt.grid()
	plt.plot(v_time,v_vol,'r')
	plt.xlabel('Time (t) [ms]')
	plt.ylabel('Bladder Volume (V) [ml]')

	# Plot pressure ---------------------------------
	import pandas as pd
	FR_pgn = pd.read_csv('PGN_freqs.csv')
	FR_pgn = np.array(FR_pgn)

	# Grill function for polynomial fit according to PGN firing rate
	# Grill, et al. 2016
	def fire_rate(x):
		f = 2.0E-03*x**3 - 3.3E-02*x**2 + 1.8*x - 0.5
		return f

	# Grill function for polynomial fit according to bladder volume
	# Grill, et al. 2016
	def blad_eq(vol):
		f = 1.5*vol - 10
		return f

	# Get values for firing rate and bladder volume functions
	FR_pgn_f = fire_rate(FR_pgn)
	vol_b = blad_eq(v_vol)

	# Because we have more values for FR than volume, we need to only use
	# the FR values at the times that we have bladder volume measurements for
	FR_pgn_plotf = []
	for n in v_time:
		FR_pgn_plotf.append(FR_pgn_f[n])

	# Grill function returning pressure in units of cm H20
	# Grill, et al. 2016
	def pressure(fr,v):
		p = []
		for n in range(len(fr)):
			p_n = fr[n] + v[n]

			# Round negative pressure up to 0
			if p_n < 0:
				p_n = 0

			p.append(p_n)

		return p

	p = pressure(FR_pgn_plotf,vol_b)

	plt.figure()
	plt.grid()
	plt.plot(v_time,p,'g')
	plt.xlabel('Time (t) [ms]')
	plt.ylabel('Bladder Pressure (P) [cm H20]')

	# Plot bladder afferent firing rate -------------

	# Grill function returning bladder afferent firing rate in units of Hz
	# Grill, et al. 2016
	def blad_aff_fr(p):
		fr = []
		for n in range(len(p)):
			fr_n = -3.0E-08*p[n]**5 + 1.0E-5*p[n]**4 - 1.5E-03*p[n]**3 + 7.9E-02*p[n]**2 - 0.6*p[n]
			
			# Round negative firing rate up to 0
			if fr_n < 0:
				fr_n = 0

			fr.append(fr_n)

		return fr 

	bladaff_fr = blad_aff_fr(p)

	plt.figure() 
	plt.grid()
	plt.plot(v_time,bladaff_fr,'m')
	plt.xlabel('Time (t) [ms]')
	plt.ylabel('Bladder Afferent Firing Rate [Hz]')
	plt.show()
