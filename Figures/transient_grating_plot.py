import numpy as np
import matplotlib.pyplot as plt
import matplotlib
#from tqdm import tqdm
from scipy.ndimage import uniform_filter1d
from scipy.interpolate import interp1d
#from scipy.optimize import curve_fit
#from scipy import signal
matplotlib.use('Qt5Agg')
plt.ion()

def interpolation(x, y, new_x):
    interpolator = interp1d(x, y, kind="linear")
    return interpolator(new_x)

pathname = '/Users/jeffreysuen/Photoacoustics Research/Data/transient grating measurement/Run00007_delta_data.txt'
delta_data = np.genfromtxt(pathname)
delays = delta_data[:,0]
delta_averaged = delta_data[:,1]
delta_0 = delta_data[:,2]
delta_1 = delta_data[:,3]
delta_2 = delta_data[:,4]

delays_even = np.arange(delays[0],delays[-1],1)
delta_0 = interpolation(delays, delta_0, delays_even)
delta_1 = interpolation(delays, delta_1, delays_even)
delta_2 = interpolation(delays, delta_2, delays_even)

halfwidth = 3.1
halfheight = 2.4

plt.figure(1, figsize=(halfwidth,3))
plt.plot(delays_even, delta_0)
plt.plot(delays_even, delta_1)
plt.plot(delays_even, delta_2)
plt.xlabel('delays[ps]')
plt.ylabel(r'$\Delta\eta/\eta$')
plt.legend(['scan 1','scan 2','scan 3'])
plt.tight_layout()
plt.show()

input('Enter to exit...')
plt.close()