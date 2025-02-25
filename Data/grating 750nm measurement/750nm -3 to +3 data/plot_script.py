import numpy as np
import matplotlib.pyplot as plt
import matplotlib as matpl
from tqdm import tqdm
from scipy.ndimage import uniform_filter1d
from scipy.optimize import curve_fit
from scipy import signal
import os

abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
#dname = os.path.abspath('') # for ipynb
os.chdir(dname)

plotlist = [1,2,3,4,5,6,7]
Abs = 0
Exp_decay = 0
Savedata = 0

peak_index = 49
delays_uncut = np.genfromtxt('Run00004/result.txt', skip_header=67)[:,1]
delta_list = []
delta_sd_list = []
run_list = [10, 7, 6, 4, 5, 8, 9]
angle_label = ['+3','+2','+1','0','-1','-2','-3']

# Read data, calculate delta and write into file
def readWriteDelta():
    N_samples = 2000
    iter_angle = iter(angle_label)
    for run in run_list:
        raw_data = np.genfromtxt(f'Run000{run:02}.txt')
        signal_data = raw_data[0,:]
        reference_data = raw_data[1,:]
        
        signal_pumped = np.reshape(signal_data[::2], (-1, N_samples//2))
        signal_blocked = np.reshape(signal_data[1::2], (-1, N_samples//2))
        reference_pumped = np.reshape(reference_data[::2], (-1, N_samples//2))
        reference_blocked = np.reshape(reference_data[1::2], (-1, N_samples//2))
        signal_pumped_averaged = np.mean(signal_pumped, axis=1)
        signal_blocked_averaged = np.mean(signal_blocked, axis=1)
        reference_pumped_averaged = np.mean(reference_pumped, axis=1)
        reference_blocked_averaged = np.mean(reference_blocked, axis=1)

        # absolute difference
        if Abs == 1: 
            delta = signal_data[::2] - signal_data[1::2] * (reference_data[::2]/reference_data[1::2])
            difflabel = 'Absolute [V]'
        # relative difference
        elif Abs == 0: 
            delta = signal_data[::2]/signal_data[1::2] * reference_data[1::2]/reference_data[::2] - 1
            difflabel = 'Relative'

        delta = np.reshape(delta, (-1, N_samples//2))
        delta_averaged = np.mean(delta, axis=1)
        delta_sd = np.std(delta, axis=1)
        delta_list.append(delta_averaged)
        delta_sd_list.append(delta_sd)

        plt.figure()
        plt.plot(signal_pumped_averaged)
        plt.plot(signal_blocked_averaged)
        plt.plot(reference_pumped_averaged)
        plt.plot(reference_blocked_averaged)
        plt.ylim(6.45,6.65)
        plt.legend(['signal pumped', 'signal blocked', 'reference pumped', 'reference blocked'])
        plt.ylabel('V')
        plt.title(f'Raw data at {next(iter_angle)}')
        plt.show()
        print(f'Read {run} done')
    np.savetxt('delta_various_angles.txt', np.vstack(delta_list))
    np.savetxt('delta_sd_various_angles.txt', np.vstack(delta_sd_list))

# readWriteDelta()

delays = delays_uncut[peak_index:]
delta_arr_uncut = np.genfromtxt('delta_various_angles.txt')
delta_arr_uncut[4,:] += 0.006
delta_sd = np.genfromtxt('delta_sd_various_angles.txt')
# mean removed
delta_arr = delta_arr_uncut[:,peak_index:]
delta_meanrm = delta_arr - np.mean(delta_arr, axis=1)[:, np.newaxis]
# FFT delta
dt = 2e-12
cutoff_freq_inGHz = 4 
cutoff_freq = cutoff_freq_inGHz * 1e9
normalized_cutoff = 2 * cutoff_freq * dt
sos = signal.butter(6, normalized_cutoff, 'hp', output='sos')
filtered_delta = np.array([signal.sosfilt(sos, row) for row in delta_meanrm])
#filtered_delta = delta_sa_bgremoved
outputlen = 2**15
delta_fft_filtered = np.abs(np.fft.rfft(filtered_delta, axis=1, n=outputlen))
delta_fft_unfiltered = np.abs(np.fft.rfft(delta_meanrm, axis=1, n=outputlen))
delta_freq = np.fft.rfftfreq(outputlen, d=dt) * 1e-9  # in Ghz, shifted starting with 0

# Standard deviation
if 7 in plotlist:
    plt.figure(7)
    color = iter(plt.cm.rainbow(np.linspace(0, 1, 7)))
    for i in range(7):
        plt.plot(delays_uncut, delta_sd[i,:], label=fr'${angle_label[i]}\degree$', c=next(color))
    plt.xlabel('Delay (ps)')
    plt.ylabel('Standard deviation (V)')
    plt.legend(loc='upper right')
    plt.title(fr'Delta trace at different angles')
    plt.show(block=False)

# Original traces
if 1 in plotlist:
    plt.figure(1)
    color = iter(plt.cm.rainbow(np.linspace(0, 1, 7)))
    for i in range(7):
        plt.plot(delays_uncut, delta_arr_uncut[i,:], label=fr'${angle_label[i]}\degree$', c=next(color))
    plt.xlabel('Delay (ps)')
    plt.ylabel('Relative')
    plt.legend(loc='upper right')
    plt.title(fr'Delta trace at different angles')
    plt.show(block=False)

# Mean removed
if 2 in plotlist:
    plt.figure(2)
    color = iter(plt.cm.rainbow(np.linspace(0, 1, 7)))
    for i in range(7):
        plt.plot(delays, delta_meanrm[i,:], label=fr'${angle_label[i]}\degree$', c=next(color))
    plt.xlabel('Delay (ps)')
    plt.ylabel('Relative')
    plt.legend(loc='upper right')
    plt.title(fr'Delta trace at different angles, mean removed')
    plt.show(block=False)

# Delta unfiltered
if 3 in plotlist:
    plt.figure(3)
    color = iter(plt.cm.rainbow(np.linspace(0, 1, 7)))
    filtered_delta_smooth = uniform_filter1d(filtered_delta, axis=1, size=3)
    max_mag = np.max(filtered_delta)*1.2
    offset = np.arange(0, 7 * max_mag, max_mag)
    for i in range(7):
        plt.plot(delays, delta_meanrm[i,:] + offset[i], label=fr'${angle_label[i]}\degree$', c=next(color))
    plt.xlabel('Delay [ps]')
    plt.ylabel('Relative')
    # put on grid lines for each trace
    plt.yticks(offset)
    plt.grid(axis='y')
    # place legend outside figure
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1)) 
    plt.subplots_adjust(right=0.85)
    plt.title(fr'Delta trace at different angles, unfiltered')
    plt.show(block=False)

# Delta filtered
if 4 in plotlist:
    plt.figure(4)
    color = iter(plt.cm.rainbow(np.linspace(0, 1, 7)))
    filtered_delta_smooth = uniform_filter1d(filtered_delta, axis=1, size=3)
    max_mag = np.max(filtered_delta)*1.2
    offset = np.arange(0, 7 * max_mag, max_mag)
    for i in range(7):
        plt.plot(delays, filtered_delta[i,:] + offset[i], label=fr'${angle_label[i]}\degree$', c=next(color))
    plt.xlabel('Delay [ps]')
    plt.ylabel('Relative')
    # put on grid lines for each trace
    plt.yticks(offset)
    plt.grid(axis='y')
    # place legend outside figure
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1)) 
    plt.subplots_adjust(right=0.85)
    plt.title(fr'Delta trace at different angles, {cutoff_freq_inGHz} Ghz HP filter')
    plt.show(block=False)
    
# FFT filtered
if 5 in plotlist:
    plt.figure(5)
    color = iter(plt.cm.rainbow(np.linspace(0, 1, 7)))
    max_mag = np.max(delta_fft_filtered)*1.05
    offset = np.arange(0, 7 * max_mag, max_mag)
    for i in range(7):
        plt.plot(delta_freq, delta_fft_filtered[i,:] + offset[i], label=fr'${angle_label[i]}\degree$', c=next(color))
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('Amplitude (Ab. Unit)')
    plt.xlim(0,60)
    plt.yticks(offset)
    plt.grid(axis='y')
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
    plt.subplots_adjust(right=0.85)
    plt.title(fr'FFT frequencies, {cutoff_freq_inGHz} Ghz HP filter')
    plt.show(block=False)

# FFT unfiltered
if 6 in plotlist:
    plt.figure(6)
    color = iter(plt.cm.rainbow(np.linspace(0, 1, 7)))
    max_mag = np.max(delta_fft_unfiltered)*1.05
    offset = np.arange(0, 7 * max_mag, max_mag)
    for i in range(7):
        plt.plot(delta_freq, delta_fft_unfiltered[i,:] + offset[i], label=fr'${angle_label[i]}\degree$', c=next(color))
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('Amplitude (Ab. Unit)')
    plt.xlim(0,60)
    plt.yticks(offset)
    plt.grid(axis='y')
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
    plt.subplots_adjust(right=0.85)
    plt.title(fr'FFT frequencies, unfiltered')
    plt.show(block=False)

input("Press Enter to continue...")



print(' ')