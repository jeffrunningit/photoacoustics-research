import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib
import glob
#from tqdm import tqdm
from scipy.ndimage import uniform_filter1d
#from scipy.interpolate import interp1d
#from scipy.optimize import curve_fit
#from scipy import signal
matplotlib.use('Qt5Agg')
plt.ion()

directory = '/Users/jeffreysuen/Photoacoustics Research/Data/grating spectrum/18012024'
s_flat_filenames = sorted(glob.glob(f'{directory}/s*flat*'))
s_grat_filenames = sorted(glob.glob(f'{directory}/s*grat*'))
p_flat_filenames = sorted(glob.glob(f'{directory}/p*flat*'))
p_grat_filenames = sorted(glob.glob(f'{directory}/p*grat*'))
for name in s_flat_filenames:
    print(name)

angles = np.arange(31,44)
legend = np.char.add(angles.astype(str), r'$\degree$')

wavelengths = np.loadtxt(s_flat_filenames[0], skiprows=18)[:-4,0]

s_flat = []
for filename in s_flat_filenames:
    s_flat.append(np.loadtxt(filename, skiprows=18)[:-4,1])
s_flat = np.array(s_flat)

s_grat = []
for filename in s_grat_filenames:
    s_grat.append(np.loadtxt(filename, skiprows=18)[:-4,1])
s_grat = np.array(s_grat)

p_flat = []
for filename in p_flat_filenames:
    p_flat.append(np.loadtxt(filename, skiprows=18)[:-4,1])
p_flat = np.array(p_flat)

p_grat = []
for filename in p_grat_filenames:
    p_grat.append(np.loadtxt(filename, skiprows=18)[:-4,1])
p_grat = np.array(p_grat)

p_gratflat = p_grat/p_flat
p_gratflat = uniform_filter1d(p_gratflat, size=3)

ps_grat = p_grat/s_grat
ps_grat = uniform_filter1d(ps_grat, size=3)

plt.figure(1)
c = iter(cm.rainbow(np.linspace(0,1,13)))
angle = 31
for spectrum in p_gratflat:
    if angle == 37:
        plt.plot(wavelengths, spectrum, color='black')
    else:
        plt.plot(wavelengths, spectrum, color=next(c))
    angle += 1
plt.xlabel(r'wavelength [nm]')
plt.ylabel(r'ratio')
plt.xlim(350,980)
plt.ylim(0,0.8)
plt.title('p grat/flat')
plt.legend(legend)
#plt.legend(['scan 1','scan 2','scan 3'])
plt.tight_layout()
plt.show()

plt.figure(2)
c = iter(cm.rainbow(np.linspace(0,1,13)))
for spectrum in ps_grat:
    plt.plot(wavelengths, spectrum, color=next(c))
plt.xlabel(r'wavelength [nm]')
plt.ylabel(r'ratio')
plt.xlim(350,980)
plt.ylim(0,0.8)
plt.title('p/s ratio')
plt.legend(legend)
plt.tight_layout()
plt.show()

plt.figure(3)
plt.plot(angles, ps_grat[:,np.argmin(abs(wavelengths - 750))])
plt.xlabel(r'Angle of incidence [deg]')
plt.ylabel(r'ratio')
#plt.xlim(350,980)
#plt.ylim(0,0.8)
plt.title('p/s ratio')
#plt.legend(legend)
plt.tight_layout()
plt.show()

# plot spectrum at 37 deg after pumped
data = np.load('/Users/jeffreysuen/Photoacoustics Research/Data/grating whitelight measurement/Measurements/data.npz')
delays = data['delays']
wavelengths_pumped = data['wavelengths']
sig = data['signal']# (pixels, delays, 4 runs)
sig_mean = sig.mean(2)[:,100] # in percentage

spectrum37 = ps_grat[6]
spectrum37 = np.interp(wavelengths_pumped, wavelengths, spectrum37)
spectrum37_pumped = spectrum37 * (1+sig_mean * 20)
plt.figure(4)
plt.plot(wavelengths_pumped, spectrum37)
plt.plot(wavelengths_pumped, spectrum37_pumped)
#plt.plot(wavelengths_pumped, sig_mean)
plt.xlabel('Wavelength (nm)')
plt.ylabel('Reflectance ratio')
plt.legend(['Before pump', 'After pump'])
plt.title('Change exaggerated 20x')
plt.tight_layout()
plt.show()

plt.figure(5)
plt.plot(wavelengths_pumped, spectrum37)
plt.plot(wavelengths_pumped, spectrum37 * (1+sig_mean))
plt.plot(wavelengths_pumped, sig_mean*100)
plt.xlabel('Wavelength (nm)')
plt.ylabel('Reflectance ratio')
plt.legend(['Reflection spectrum unpumped', 'Reflection spectrum pumped', 'Pump-induced change'])
plt.title('Change exaggerated 100x')
plt.tight_layout()
plt.show()

from scipy.interpolate import CubicSpline

# Fit a cubic spline to the data
spline = CubicSpline(wavelengths_pumped[::-1], spectrum37[::-1])

# Generate points for plotting
x_fit = np.linspace(wavelengths_pumped[-1], wavelengths_pumped[0], 100)
y_fit = spline(x_fit)
y_derivative = spline.derivative()(x_fit)

# Plot the original data, spline, and its derivative
plt.figure(figsize=(10, 6))
plt.scatter(wavelengths_pumped, spectrum37, label='Data', color='red', marker='.')
plt.plot(x_fit, y_fit, label='Cubic Spline Fit', color='blue')
plt.plot(x_fit, y_derivative, label='Derivative of Cubic Spline', linestyle='--', color='green')
plt.plot(wavelengths_pumped, sig_mean)
plt.xlabel('wavelength (nm)')
plt.ylabel('y')
plt.title('Cubic Spline Fitting and Derivative Calculation')
plt.legend()
plt.show()

input('Enter to exit...')
print('.')