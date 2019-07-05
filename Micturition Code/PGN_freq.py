import pandas as pd 
import numpy as np

PGN_data = pd.read_csv('data.csv')
current_data = []
#print(PGN_data)

low = 0
high = 1000
freqs = []
tf = []

for i in np.arange(1,10):

    for n in range(len(PGN_data)):
        if PGN_data[n] >= low and PGN_data[n] <= high:
            tf.append(True)
        else:
            tf.append(False)

    current_data = PGN_data[tf]
    freq = len(current_data)/1000
    freqs.append(freq)
        
    low += 1000
    high += 1000

print(freqs)