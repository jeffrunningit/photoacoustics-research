import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
from scipy.ndimage import uniform_filter1d
from scipy.optimize import curve_fit
from scipy import signal

#matplotlib.use("MacOSX")

os.chdir("7 angles data")
allfilenames = os.listdir()

peak_index = 40 # Starting peak position (2ps) where mean is calculated from

Run_indices =  [20, #p1
                18, #p0.5
                22, #p0.25
                16, #0
                21, #m0.25
                17, #m0.5
                19] #m1
# in order of angles
plotname = 'delta_various_angles_smooth'
plotlist = [1,2,3,4,5,6,7]
Savefig = 0

Run_indices = [str(n).zfill(5) for n in Run_indices]
plotfilenames = []
for n in Run_indices:
    for file in allfilenames:
        if n in file:   
            plotfilenames.append(file)
plotdata = [np.genfromtxt(file) for file in plotfilenames]

delays = plotdata[0][:,0]
delta_interp_list = []
background_list = []
delta_smooth_list = []
sin_popt = []
delta_sinremove = []
delta_freqfiltered = []
fft_list = []
freqs_list = []
filter_size = 3
timestep = 0.1
cutoff_freq = 4e9
sos = signal.butter(4, cutoff_freq, 'hp', fs=1/0.1e-12, output='sos')

def sinwave(t, a, p, T, c):
    return a * np.sin(2 * np.pi / p * t + T) + c

for filename, data in zip(plotfilenames, plotdata):
    delta = data[:,1] * 100 # change to percentage %

    # interpolation
    delays_uniform = np.arange(data[peak_index,0], max(data[:,0]),timestep)
    delta_interp = np.interp(delays_uniform, delays, delta)
    delta_interp_list.append(delta_interp)

    # calculate background mean for each run (interpolated)
    delta_mean = np.mean(delta_interp)
    background_list.append(delta_mean)

    # apply smooth filter
    delta_smooth = uniform_filter1d(delta, size=filter_size)
    delta_smooth_list.append(delta_smooth)

    # fit to cosine curve
    popt, pcov = curve_fit(sinwave, delays, delta, p0=[max(delta)-delta_mean, 350, 0, delta_mean])
    sin_popt.append(popt)
    delta_sinremove.append(delta - sinwave(delays, *popt))

    # highpass filter
    filtered = signal.sosfilt(sos, delta_interp[peak_index:])
    delta_freqfiltered.append(filtered)

    # FFT
    outputlen = 2**15
    delta_fft = np.fft.fftshift(np.fft.fft(filtered, n=outputlen))
    delta_freqs = np.fft.fftshift(np.fft.fftfreq(outputlen, d=0.1e-12) * 1e-9)  # in Ghz, shifted starting with negative
    fft_list.append(delta_fft)
    freqs_list.append(delta_freqs)



figsize = (14,7)


# 1. Simple plot 
if 1 in plotlist:
    plt.figure(figsize=figsize)
    color = iter(plt.cm.rainbow(np.linspace(0, 1, len(plotfilenames))))
    for filename, delta in zip(plotfilenames, delta_smooth_list):
        plt.plot(delays, delta, label=filename, c=next(color))
    plt.xlabel('Delay (ps)')
    plt.ylabel(r'$\Delta R/R_0$ (%)')
    plt.title('Signal from multiple runs')
    plt.legend()
    plt.show(block=False)

# 2. Mean removed
if 2 in plotlist:
    plt.figure(figsize=figsize)
    color = iter(plt.cm.rainbow(np.linspace(0, 1, len(plotfilenames))))
    for filename, delta, background in zip(plotfilenames, delta_smooth_list, background_list):
        plt.plot(delays, delta-background, label=filename, c=next(color))
    plt.xlabel('Delay (ps)')
    plt.ylabel(r'$\Delta R/R_0$ (%)')
    plt.title('Signal from multiple runs (background mean removed)')
    plt.legend()
    plt.show(block=False)

# 3. Offset
if 3 in plotlist:
    plt.figure(figsize=figsize)
    offset_list = np.arange(0,len(Run_indices) * 0.2, 0.2)
    offset = iter(offset_list)
    color = iter(plt.cm.rainbow(np.linspace(0, 1, len(plotfilenames))))
    for filename, delta, background, popt in zip(plotfilenames, delta_smooth_list, background_list, sin_popt):
        c = next(color)
        off = next(offset)
        plt.plot(delays, delta - background + off, label=filename, c=c)
        plt.plot(delays, sinwave(delays, *popt) - background + off, linestyle=':', c=c)
    plt.yticks(offset_list)
    plt.grid(axis='y')
    plt.xlabel('Delay (ps)')
    plt.ylabel(r'$\Delta R/R_0$ (%)')
    plt.title('Signal from multiple runs (offset centered at background mean)')
    plt.legend()
    plt.show(block=False)

# 4. Sine removed
if 4 in plotlist:
    plt.figure(figsize=figsize)
    offset_list = np.arange(0,len(Run_indices) * 0.2, 0.2)
    offset = iter(offset_list)
    color = iter(plt.cm.rainbow(np.linspace(0, 1, len(plotfilenames))))
    for filename, delta in zip(plotfilenames, delta_sinremove):
        plt.plot(delays, delta + next(offset), label=filename, c=next(color))
    plt.yticks(offset_list)
    plt.grid(axis='y')
    plt.xlabel('Delay (ps)')
    plt.ylabel(r'$\Delta R/R_0$ (%)')
    plt.title('Signal from multiple runs (sin wave removed)')
    plt.legend()
    plt.show(block=False)

# 5. frequency filter highpass offset
if 5 in plotlist:
    plt.figure(figsize=figsize)
    offset_list = np.arange(0,len(Run_indices) * 0.2, 0.2)
    offset = iter(offset_list)
    color = iter(plt.cm.rainbow(np.linspace(0, 1, len(plotfilenames))))
    for filename, delta in zip(plotfilenames, delta_freqfiltered):
        plt.plot(delays_uniform[peak_index:], delta + next(offset), label=filename, c=next(color))
    plt.yticks(offset_list)
    plt.grid(axis='y')
    plt.xlabel('Delay (ps)')
    plt.ylabel(r'$\Delta R/R_0$ (%)')
    plt.title(f'Signal from multiple runs ({cutoff_freq/1e9:.2g} GHz HP)')
    plt.legend()
    plt.show(block=False)

# 6. frequency space
if 6 in plotlist:
    plt.figure(figsize=figsize)
    color = iter(plt.cm.rainbow(np.linspace(0, 1, len(plotfilenames))))
    offset = 0
    offset_spacing = np.max(np.abs(fft_list[3]))*1.1
    offset_list = np.arange(0, len(plotfilenames) * offset_spacing, offset_spacing)
    for filename, freq, fft in zip(plotfilenames, freqs_list, fft_list):
        plt.plot(freq, np.abs(fft) + offset, label=filename, c=next(color))
        offset += offset_spacing 
    plt.yticks(offset_list)
    plt.grid(axis='y')
    plt.xlim(0,60)
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('Amplitude (Arb. Unit)')
    plt.title(f'Frequencies from multiple runs ({cutoff_freq/1e9:.2g} GHz HP)')
    #plt.legend(loc='upper right')
    plt.show(block=False)

# 7. frequency filter highpass
if 7 in plotlist:
    plt.figure(figsize=figsize)
    color = iter(plt.cm.rainbow(np.linspace(0, 1, len(plotfilenames))))
    for filename, delta in zip(plotfilenames, delta_freqfiltered):
        plt.plot(delays_uniform[peak_index:], delta, label=filename, c=next(color))
    plt.grid(axis='y')
    plt.xlabel('Delay (ps)')
    plt.ylabel(r'$\Delta R/R_0$ (%)')
    plt.title(f'Signal from multiple runs ({cutoff_freq/1e9:.2g} GHz HP)')
    plt.legend()
    plt.show(block=False)

if Savefig == 1:
    print('saving figures')
    os.chdir('..\combined_plots')
    if 1 in plotlist:
        fig = plt.figure(plotlist.index(1)+1)
        plt.savefig(plotname + f'.png', dpi=200)
    if 2 in plotlist:
        fig = plt.figure(plotlist.index(2)+1)
        plt.savefig(plotname + f'_bgremoved.png', dpi=200)
    if 3 in plotlist:
        fig = plt.figure(plotlist.index(3)+1)
        plt.savefig(plotname + f'_offset.png', dpi=200)
    if 4 in plotlist:
        fig = plt.figure(plotlist.index(4)+1)
        plt.savefig(plotname + f'_sineremoved.png', dpi=200)
    if 5 in plotlist:
        fig = plt.figure(plotlist.index(5)+1)
        plt.savefig(plotname + f'_hp_time_offset.png', dpi=200)
    if 6 in plotlist:
        fig = plt.figure(plotlist.index(6)+1)
        plt.savefig(plotname + f'_hp_freq.png', dpi=200)
    if 7 in plotlist:
        fig = plt.figure(plotlist.index(7)+1)
        plt.savefig(plotname + f'_hp_time.png', dpi=200)
    print('saved figures')



input('Enter to exit...')