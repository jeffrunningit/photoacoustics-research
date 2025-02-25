import numpy as np
import matplotlib.pyplot as plt
import oberon_package.dataprocessing as ob
import oberon_package.functions as fn
import matplotlib as matpl
from oberon_package.parameters import Options
from tqdm import tqdm
from scipy.ndimage import uniform_filter1d
from scipy.optimize import curve_fit
from scipy import signal

import os

Run = "Run00010"
angle = '0'    # relative to resonance angle. p=plus, m=minus, 0=resonance
Abs = 0         # 1=absolute difference, 0=relative difference
Exp_decay = 0   # 1=remove exponential decay, 0=remove average value
peak_index = 40       # index of first peak (t=2ps)
Savedata = 0    # Save data in files
Savefig = 0    # Save figures as png
plotlist = [1,3,4,5,6]


read_data = ob.ReadPlus(folder = Run,dname=dname)
data = read_data.read_digitizer_data()

raw_data = data['raw_data']
signal_data = raw_data['Signal']
reference_data = raw_data['Reference']
delays = np.array(data['delays'])
delays_order = data['delays_order']
N_scans = data['scans']
N_samples = data['samples']

raw_data_toSave = np.vstack([signal_data, reference_data])
np.savetxt(Run+'.txt', raw_data_toSave)

# Even data = pumped, Odd data = blocked

signal_pumped = np.reshape(signal_data[::2], (-1, N_samples//2))
signal_blocked = np.reshape(signal_data[1::2], (-1, N_samples//2))
reference_pumped = np.reshape(reference_data[::2], (-1, N_samples//2))
reference_blocked = np.reshape(reference_data[1::2], (-1, N_samples//2))
averaged_pumped_signal = np.mean(signal_pumped, axis=1)
averaged_blocked_signal = np.mean(signal_blocked, axis=1)
averaged_pumped_ref = np.mean(reference_pumped, axis=1)
averaged_blocked_ref = np.mean(reference_blocked, axis=1)


# absolute difference
if Abs == 1: 
    delta = signal_data[::2] - signal_data[1::2] * (reference_data[::2]/reference_data[1::2])
    difflabel = 'Absolute [V]'
# relative difference
elif Abs == 0: 
    delta = signal_data[::2]/signal_data[1::2] * reference_data[1::2]/reference_data[::2] - 1
    difflabel = 'Relative'


delta = np.reshape(delta, (-1, N_samples//2))
delta_shots_averaged = np.mean(delta, axis=1)
delta_scans = np.reshape(delta_shots_averaged, (N_scans,-1))
delays_indices = np.argsort(np.reshape(delays_order, (N_scans,-1)))
for i in range(N_scans):
    delta_scans[i] = delta_scans[i][delays_indices[i]]
delta_scans_averaged = np.mean(delta_scans, axis=0)


# interpolation
timestep = 0.1
delays_uniform = np.arange(delays[peak_index],max(delays),timestep)
delta_scans_averaged_interp = np.interp(delays_uniform, delays[peak_index:], delta_scans_averaged[peak_index:])

# subtract exponential decay
if Exp_decay == 1:
    def exp_func(x, a, b, c):
        return a * np.exp(-b * x) + c
    param = (-0.0021, 0.007, 0)
    param, cov = curve_fit(exp_func, delays[peak_index:], delta_scans_averaged[peak_index:], param)
    delta_bg = exp_func(delays_uniform, *param)
    delta_sa_bgremoved = delta_scans_averaged_interp - delta_bg
# subtract mean
elif Exp_decay == 0:
    delta_sa_bgremoved = delta_scans_averaged_interp - np.mean(delta_scans_averaged_interp)


# FFT delta
d = timestep * 1e-12
sos = signal.butter(10, 4e9, 'hp', fs=1/d, output='sos')
filtered_delta = signal.sosfilt(sos, delta_scans_averaged_interp)
#filtered_delta = delta_sa_bgremoved
outputlen = 2**17
delta_fft = np.fft.fftshift(np.fft.fft(filtered_delta, n=outputlen))
delta_freq = np.fft.fftshift(np.fft.fftfreq(outputlen, d=d)*1e-9)  # in Ghz, shifted starting with 0


# Saving data into txt file for later
if Savedata == 1:
    avg_data = np.transpose(np.array([delays_order, averaged_pumped_signal, averaged_blocked_signal, averaged_pumped_ref, averaged_blocked_ref]))
    header = 'delays[ps] pumped_signal, blocked_signal, pumped_reference, blocked_reference'
    np.savetxt(f'{Run}_{angle}deg_avg_data_overshots.txt', avg_data, header=header)

    delta_data = np.column_stack((delays, delta_scans_averaged, delta_scans.T))
    header = 'delays[ps] delta_averaged'
    for i in range(N_scans):
        header += f' scan {i}'
    np.savetxt(f'{Run}_{angle}deg_delta_data.txt', delta_data, header=header)

# Plotting all pumped and blocked signal and reference in sequence of measurement
# averaged every time step
if 1 in plotlist:
    plt.figure(1)
    plt.plot(averaged_blocked_signal, label = 'Blocked signal')
    plt.plot(averaged_pumped_signal, label='Pumped signal')
    plt.plot(averaged_blocked_ref, label = 'Blocked reference')
    plt.plot(averaged_pumped_ref, label='Pumped reference')
    plt.legend()
    plt.xlabel('Data points')
    plt.ylabel(r'V')
    plt.title(f"{Run} at {angle} deg all averaged data")
    plt.show(block=False)

# Plotting every run delta with respect to delay time
# averaged every time step
smoothsize = 3
if 2 in plotlist:
    smooth_delta_scans = np.array([uniform_filter1d(data, size=smoothsize) for data in delta_scans])
    plt.figure(2)
    for i in range(N_scans):
        plt.plot(delays, smooth_delta_scans[i], label = 'scan '+str(i+1))
    plt.legend()
    plt.xlabel('Delay (ps)')
    plt.ylabel(difflabel)
    plt.title(fr'{Run} at {angle} deg $\eta$ (smooth {smoothsize})')
    plt.show(block=False)

# Plotting delta with respect to delay time
# averaged between runs
if 3 in plotlist:
    smooth_delta_scans_averaged = uniform_filter1d(delta_scans_averaged, size=smoothsize)
    plt.figure(3)
    plt.plot(delays, smooth_delta_scans_averaged)
    if Exp_decay == 1:
        plt.plot(delays_uniform, delta_bg)
    plt.legend("Fitted exp decay", "Data")
    plt.xlabel('Delay (ps)')
    plt.ylabel(difflabel)
    plt.title(fr'{Run} at {angle} deg $\eta$ {N_scans} runs averaged (smooth {smoothsize})')
    plt.show(block=False)

# Plotting delta by data points
# averaged between runs
if 4 in plotlist:
    plt.figure(4)
    plt.plot(delta_scans_averaged, label = 'signal')
    plt.legend()
    plt.xlabel('Data points')
    plt.ylabel(difflabel)
    plt.title(fr'{Run} at {angle} deg $\eta$ {N_scans} runs averaged')
    plt.show(block=False)

# FFT of delta
if 5 in plotlist:
    plt.figure(5)
    plt.plot(delta_freq, np.abs(delta_fft))
    plt.xlim(0,20)
    plt.legend()
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('Amplitude (Arb. Unit)')
    plt.title(fr'{Run} at {angle} deg $\eta$ FFT')
    plt.show(block=False)

# Thermal background removed
if 6 in plotlist:
    plt.figure(6)
    plt.plot(delays_uniform, delta_sa_bgremoved, label='Delta')
    plt.plot(delays_uniform, filtered_delta, label='Frequency filtered')
    plt.legend()
    plt.xlabel('Delay (ps)')
    plt.ylabel(difflabel)
    plt.title(fr'{Run} at {angle} deg $\eta$ background removed')
    plt.show(block=False) 

if Savefig == 1:
    print('saving figures')
    for i in plt.get_fignums():
        fig = plt.figure(i)
        fig.set_size_inches(16, 9)
        plt.savefig(f'{Run}_{angle}deg_fig{i}.png', dpi=200)
    print('saved figures')
print('.')



print(' ')