import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import glob

'''
from inspect import getsourcefile
from os.path import abspath
print(f'The directory of current file: {abspath(getsourcefile(lambda:0))}')
'''

directory = '/Users/jeffreysuen/Photoacoustics Research/Data/grating 750nm measurement/Run00001 1000ps'
filenames = np.sort(glob.glob(f'{directory}/*delta_data*'))

delays_full = np.loadtxt(filenames[0])[:,0]
data = []
for filename in filenames:
    run_data = np.loadtxt(filename)[:,1]
    run_data -= run_data[5]
    data.append(run_data)
data = np.vstack(data) * 100
# Saving data into npz file for later
np.savez(f'{directory}/processed_data', delays_full=delays_full, data=data, cutindex=49)

print('.')