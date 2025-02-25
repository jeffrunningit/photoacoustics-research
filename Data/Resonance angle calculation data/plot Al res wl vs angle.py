import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy.interpolate as interp

matplotlib.use('Qt5Agg')
plt.ion()

abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)


filename = "/Users/jeffreysuen/Photoacoustics Research/Data/Resonance angle calculation data/Al_diffraction_angles.csv"
data = np.genfromtxt(filename, delimiter=',')[91:152,:]
wl = data[:,0]
n = data[:,1]
k = data[:,2]
theta = -data[:,3]

plt.figure(1).set_size_inches(4.6,4)
plt.plot(theta, wl*1e9)
plt.ylim(650,850)
plt.xlim(28,38)
major_xticks = np.arange(28,39,2)
minor_xticks = np.arange(28,38)
major_yticks = np.arange(650,851,50)
minor_yticks = np.arange(650,850,25)
plt.xticks(major_xticks, [f'{tick:.0f}°' for tick in major_xticks])
plt.xticks(minor_xticks, minor=True)
plt.yticks(major_yticks)
plt.yticks(minor_yticks, minor=True)
plt.grid(which='major', alpha=0.5)
plt.grid(axis='x', which='minor', alpha=0.2)
plt.ylabel(r'SPP resonance wavelength (nm)')
plt.xlabel('Angle of incidence')
plt.axhline(750, color='tab:red',linewidth=1)
plt.axvline(33.1, color='tab:grey',linewidth=1)
#plt.title('Resonance wavelength for each angle')
plt.tight_layout()
plt.show()

c = 3e8
em = n**2 - k**2
kspp = 2 * np.pi / wl * np.sqrt(em/(em+1))
kx = 2 * np.pi / wl * np.sin(38 * 2*np.pi/360)
kg = 2 * np.pi / 1600e-9

plt.figure(2)
plt.plot(kspp, 2 * np.pi * c / wl)
plt.plot(kx + kg, 2 * np.pi * c / wl)
plt.legend(['kspp', 'kx+kg'])
plt.xlabel('wavelength (nm)')
plt.ylabel(r'k vector ($m^{-1}$)')
plt.show()

x = np.linspace(-3,3,100)
y1 = 1-np.exp(-x**2)
y2 = 1-np.exp(-(x-0.7)**2)
y3 = 1-np.exp(-(x+0.7)**2)
plt.figure(3)
plt.plot(x,y1,color='black')
plt.plot(x,y2,color='blue')
plt.plot(x,y3,color='red')
plt.xticks([])
plt.yticks([])
plt.title(r'$\lambda_+$ $\lambda_-$ $\lambda_{res}$')
plt.legend([r'$+\Delta\theta$ $-\Delta\theta$ $\theta_{res}$'])
plt.show()

plt.figure(4)
plt.plot(wl, em)
plt.show()


print('.')