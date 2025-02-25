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

Run = "Run00003"
Abs = 0         # 1=absolute difference, 0=relative difference
Exp_decay = 0   # 1=remove exponential decay, 0=remove average value
peak = 49       # index of first peak (t=2ps)
Savedata = 0    # Save processed data in a file
Savefig = 0     # Save figures as png
plotlist = [1,2,3,5,6]

abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
#dname = os.path.abspath('') # for ipynb
os.chdir(dname)

read_data = ob.ReadPlus(folder = Run,dname=dname)
data = read_data.read_digitizer_data()

raw_data = data['raw_data']
signal_data = raw_data['Signal']
reference_data = raw_data['Reference']
delays = np.array(data['delays'])
delays_order = data['delays_order']
N_scans = data['scans']
N_samples = data['samples']

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
delays_uniform = np.arange(delays[peak],1000,0.1)
delta_scans_averaged_interp = np.interp(delays_uniform, delays[peak:], delta_scans_averaged[peak:])

# subtract exponential decay
if Exp_decay == 1:
    def exp_func(x, a, b, c):
        return a * np.exp(-b * x) + c
    param = (-0.0021, 0.007, 0)
    param, cov = curve_fit(exp_func, delays[peak:], delta_scans_averaged[peak:], param)
    delta_bg = exp_func(delays_uniform, *param)
    delta_sa_bgremoved = delta_scans_averaged_interp - delta_bg
# subtract mean
elif Exp_decay == 0:
    delta_sa_bgremoved = delta_scans_averaged_interp - np.mean(delta_scans_averaged_interp)

# FFT delta

delta_fft = np.fft.fft(delta_sa_bgremoved)
delta_freq = np.fft.fftfreq(len(delays_uniform))

# Saving processed data into txt file
avg_data = np.transpose(np.array([delays_order, averaged_pumped_signal, averaged_blocked_signal, averaged_pumped_ref, averaged_blocked_ref]))
header = 'delays[ps] pumped_signal, blocked_signal, pumped_reference, blocked_reference'
np.savetxt(Run+'_avg_data_overshots.txt', avg_data, header=header)

delta_data = np.column_stack((delays, delta_scans_averaged, delta_scans.T))
header = 'delays[ps] delta_averaged'
for i in range(N_scans):
    header += f' scan {i}'
np.savetxt(Run+'_delta_data.txt', delta_data, header=header)

# Plotting all pumped and blocked signal and reference in sequence of measurement
# averaged every time step
if 1 in plotlist:
    plt.figure()
    plt.plot(averaged_blocked_signal, label = 'Blocked signal')
    plt.plot(averaged_pumped_signal, label='Pumped signal')
    plt.plot(averaged_blocked_ref, label = 'Blocked reference')
    plt.plot(averaged_pumped_ref, label='Pumped reference')
    plt.legend()
    plt.xlabel('Data points')
    plt.ylabel(r'V')
    plt.title(Run+" all averaged data")
    plt.show(block=False)

# Plotting delta with respect to delay time
# averaged every time step
if 2 in plotlist:
    smooth_delta_scans = np.array([uniform_filter1d(data, size=3) for data in delta_scans])
    plt.figure()
    for i in range(N_scans):
        plt.plot(delays, delta_scans[i], label = 'scan '+str(i+1))
    plt.legend()
    plt.xlabel('Delay (ps)')
    plt.ylabel(difflabel)
    plt.title(Run+r' $\eta$')
    plt.show(block=False)


# Plotting delta with respect to delay time
# averaged between runs
if 3 in plotlist:
    plt.figure()
    plt.plot(delays, delta_scans_averaged)
    if Exp_decay == 1:
        plt.plot(delays_uniform, delta_bg)
    plt.legend("Fitted exp decay", "Data")
    plt.xlabel('Delay (ps)')
    plt.ylabel(difflabel)
    plt.title(Run+fr' $\eta$ {N_scans} runs averaged')
    plt.show(block=False)

# Plotting delta by data points
# averaged between runs
if 4 in plotlist:
    plt.figure()
    plt.plot(delta_scans_averaged, label = 'signal')
    plt.legend()
    plt.xlabel('Data points')
    plt.ylabel(difflabel)
    plt.title(Run+fr' $\eta$ {N_scans} runs averaged')
    plt.show(block=False)

# FFT of delta
if 5 in plotlist:
    plt.figure()
    plt.plot(delta_freq*1e3, np.abs(delta_fft))
    plt.xlim(0,20)
    plt.legend()
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('Amplitude (Arb. Unit)')
    plt.title(Run+fr' $\eta$ {N_scans} runs averaged')
    plt.show(block=False)

# Thermal background removed
if 6 in plotlist:
    plt.figure()
    plt.plot(delays[49:], delta_sa_bgremoved)
    #plt.legend()
    plt.xlabel('Delay (ps)')
    plt.ylabel(difflabel)
    plt.title(Run+fr' $\eta$ {N_scans} runs averaged')
    plt.show(block=False) 


if Savefig == 1:
    for i in plt.get_fignums():
        plt.figure(i)
        plt.savefig(Run+f'_fig{i}.png')
print('.')



print(' ')