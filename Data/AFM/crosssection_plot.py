import numpy as np
import matplotlib.pyplot as plt
import matplotlib

#matplotlib.use('Qt5Agg')
plt.ion()

height = np.loadtxt("/Users/jeffreysuen/Photoacoustics Research/Data/AFM/grating.0_00003_height.txt")
height_av = height.mean(axis=1)
height_av *= 1e9
height_av -= np.min(height_av) + 3
x = np.linspace(0, 8500, height.shape[0])

halfwidth = 3.1
halfheight = 2.4
plt.figure(1, figsize=(3.1,3))
plt.plot(x, height_av, color='k')
plt.xlabel(r'x [$\mu m$]')
plt.ylabel(r'z [nm]')
plt.xlim(right=5200)
plt.tight_layout()
plt.show()

input("Press Enter to continue...")

print('.')