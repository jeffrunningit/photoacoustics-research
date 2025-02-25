import numpy as np
import matplotlib.pyplot as plt

spect_data = np.genfromtxt('Al flat 35.6 deg.TXT', delimiter=';', skip_header=8)
print(spect_data)

plt.plot(spect_data[:,0], spect_data[:,4])
plt.title('Al flat 35.6 deg spectrum')
plt.ylabel('counts')
plt.xlabel('wavelength [nm]')
plt.show()
print("Done.")