import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib
import glob
#from tqdm import tqdm
from scipy.ndimage import uniform_filter1d
#from scipy.interpolate import interp1d
#from scipy.optimize import curve_fit
from scipy import signal
matplotlib.use('Qt5Agg')
plt.ion()

'''
directory = '/Users/jeffreysuen/Photoacoustics Research/Data/grating 750nm measurement/20240612'
filenames = np.sort(glob.glob(f'{directory}/*delta_data*'))

delays_full = np.loadtxt(filenames[0])[:,0]
data = []
for filename in filenames:
    run_data = np.loadtxt(filename)[:,1]
    run_data -= run_data[5]
    data.append(run_data)
data = np.array(data) * 100

delays = delays_full[49:]
data_cut = data[:,49:] 
# Saving data into npz file for later
np.savez('/Users/jeffreysuen/Photoacoustics Research/Data/grating 750nm measurement/750nm measurement processed m3 to p3', delays_full=delays_full, data=data, cutindex=49)
'''
directory_choice = 3
plotlist = [3]

if directory_choice == 0: 
    directory = '/Users/jeffreysuen/Photoacoustics Research/Data/grating 750nm measurement/750nm -1 to +1 data/0 degree data'
if directory_choice == 1: 
    directory = '/Users/jeffreysuen/Photoacoustics Research/Data/grating 750nm measurement/750nm -1 to +1 data/7 angles data'
if directory_choice == 3: 
    directory = '/Users/jeffreysuen/Photoacoustics Research/Data/grating 750nm measurement/750nm -3 to +3 data'
if directory_choice == 1000: 
    directory = '/Users/jeffreysuen/Photoacoustics Research/Data/grating 750nm measurement/Run00001 1000ps'

alldata = np.load(f'{directory}/processed_data.npz', allow_pickle=True)
delays_full = alldata['delays_full']
data = alldata['data']
cutindex = alldata['cutindex']
delays = delays_full[cutindex:]
data_cut = data[:,cutindex:]

N_traces = 1
if directory_choice == 3:
    legend = [r'$-3\degree$',
                r'$-2\degree$',
                r'$-1\degree$',
                r'$\theta_{res}$',
                r'$+1\degree$',
                r'$+2\degree$',
                r'$+3\degree$']
    N_traces = 7
elif directory_choice == 1:
    legend =   [r'$-1\degree$',
                r'$-0.5\degree$',
                r'$-0.25\degree$',
                r'$\theta_{res}$',
                r'$+0.25\degree$',
                r'$+0.5\degree$',
                r'$+1\degree$']
    N_traces = 7
elif directory_choice == 0:
    legend = ['Run 1',  
            'Run 2',
            'Run 3',
            'Run 4',
            'Run 5']
    N_traces = 5


# subtract exponential decay
data_bgremoved = data_cut - np.mean(data_cut, axis=1).reshape(-1,1)

# FFT delta
dt = 2e-12
cutoff_freq_inGHz = 4
cutoff_freq = cutoff_freq_inGHz * 1e9
normalized_cutoff = 2 * cutoff_freq * dt
sos = signal.butter(4, normalized_cutoff, 'hp', output='sos')
filtered_delta = signal.sosfilt(sos, data_bgremoved, axis=1)
outputlen = 2**15
delta_fft_filtered = np.abs(np.fft.rfft(filtered_delta, axis=1, n=outputlen))

# fft unfiltered
delta_fft_unfiltered = np.abs(np.fft.rfft(filtered_delta, axis=1, n=outputlen))
delta_freq = np.fft.rfftfreq(outputlen, d=dt) * 1e-9  # in Ghz, shifted starting with 0



# smooth
smoothsize = 3
smooth_data = uniform_filter1d(data, size=smoothsize)
smooth_data_bgremoved = uniform_filter1d(data_bgremoved, size=smoothsize)

# Plot original
if 1 in plotlist:
    plt.figure(1, figsize=(6.4,4))
    c = iter(cm.rainbow(np.linspace(0, 1, N_traces)))
    label = iter(legend)
    #plt.plot(delays_full, smooth_data[0,:])
    for delta in smooth_data:
        plt.plot(delays_full, delta, color=next(c), label = next(label))
    handles, labels = plt.gca().get_legend_handles_labels()
    # Reverse the order of handles and labels
    plt.legend(handles[::-1], legend[::-1])
    yt_interval = 0.2
    ylim = plt.ylim()
    yticks = np.arange(ylim[0] - ylim[0] % yt_interval + yt_interval, ylim[1] - ylim[1] % yt_interval + yt_interval/2, yt_interval)
    plt.yticks(yticks)
    plt.grid(True, alpha=0.3)
    plt.xlabel('Delay (ps)')
    plt.ylabel(r'$\Delta R/R$ (%)')
    plt.title(fr'Original data')
    plt.tight_layout()
    plt.show()

# Mean removed
if 2 in plotlist:
    smooth_data = uniform_filter1d(data_bgremoved, size=smoothsize, axis=1)
    plt.figure(2, figsize=(6.4,4))
    color = iter(plt.cm.rainbow(np.linspace(0, 1, N_traces)))
    for delta in data_bgremoved:
        plt.plot(delays, delta, c=next(color))
    yt_interval = 0.1
    ylim = plt.ylim()
    yticks = np.arange(ylim[0] - ylim[0] % yt_interval + yt_interval, ylim[1] - ylim[1] % yt_interval + yt_interval/2, yt_interval)
    plt.yticks(yticks)
    plt.grid(True, alpha=0.3)
    plt.xlabel('Delay (ps)')
    plt.ylabel(r'$\Delta R/R$ (%)')
    plt.legend(legend)
    plt.title(fr'Delta trace at different angles, mean removed')
    plt.tight_layout()
    plt.show()

# Delta unfiltered
if 3 in plotlist:
    plt.figure(3, figsize=(3.1,4.8))
    color = iter(plt.cm.rainbow(np.linspace(0, 1, N_traces)))
    #max_mag = np.max(np.abs(data_bgremoved)) * 1.2
    offset_int = 0.4
    offsetlist = np.arange(0, 7 * offset_int, offset_int)
    for i in range(N_traces):
        plt.plot(delays, data_bgremoved[i,:] + offsetlist[i], c=next(color))
    plt.xlabel('Delay (ps)')
    plt.ylabel(r'$\Delta R/R$ (%)')
    # put on grid lines for each trace
    plt.xticks([0,100,200,300,400,500])
    plt.yticks(offsetlist)
    plt.grid(True, alpha=0.3)
    plt.ylim(-0.2,2.8)
    #plt.legend(legend_deg) 
    plt.subplots_adjust(right=0.85)
    plt.title(fr'Delta trace at different angles, unfiltered')
    plt.tight_layout()
    plt.show(block=False)

# Delta filtered
if 4 in plotlist:
    plt.figure(4)
    color = iter(plt.cm.rainbow(np.linspace(0, 1, N_traces)))
    filtered_delta_smooth = uniform_filter1d(filtered_delta, axis=1, size=3)
    max_mag = np.max(filtered_delta)*1.2
    offsetlist = np.arange(0, 7 * max_mag, max_mag)[::-1]
    for i in range(7):
        plt.plot(delays, filtered_delta[i,:] + offsetlist[i], c=next(color))
    plt.xlabel('Delay (ps)')
    plt.ylabel('Relative')
    # put on grid lines for each trace
    plt.yticks([])
    plt.grid(axis='y')
    # place legend outside figure
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1)) 
    plt.subplots_adjust(right=0.85)
    plt.title(fr'Delta trace at different angles, {cutoff_freq_inGHz} Ghz HP filter')
    plt.show(block=False)
    
# FFT filtered
if 5 in plotlist:
    plt.figure(5)
    color = iter(plt.cm.rainbow(np.linspace(0, 1, N_traces)))
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
    plt.figure(6, figsize=(3.1,4.8))
    color = iter(plt.cm.rainbow(np.linspace(0, 1, N_traces)))
    max_mag = np.max(delta_fft_unfiltered)*1.05
    offset = np.arange(0, 7 * max_mag - 0.1, max_mag)+2
    for i in range(7): 
        plt.plot(delta_freq, delta_fft_unfiltered[i,:] + offset[i], c=next(color))
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('Amplitude (Ab. Unit)')
    plt.xlim(0,40)
    plt.yticks(offset)
    plt.grid(axis='y', alpha=0.3)
    #plt.gca().set_yticklabels([])
    plt.title(fr'FFT frequencies, unfiltered')
    plt.tight_layout()
    plt.show()
    print(offset)

input("Press Enter to continue...")



print(' ')