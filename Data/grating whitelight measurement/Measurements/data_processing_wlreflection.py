# %%

import numpy as np
import matplotlib.pyplot as plt
import oberon_package.dataprocessing as ob
import oberon_package.functions as fn
import matplotlib as matpl
from oberon_package.parameters import Options
from tqdm import tqdm
import scipy.signal as ssg
import photutils.centroids as centr
import photutils.segmentation as segm

from scipy import signal
from scipy.ndimage import uniform_filter1d
from scipy import interpolate
from scipy.optimize import curve_fit

import os


options = Options()

options.processing.outlier_remove = False
options.processing.outlier_min_intensity = int(1.5e8)
options.processing.reference = "Referenced"
options.processing.calc_mode = "Relative"
options.processing.pump_order = 0
options.processing.ref_shift = 2
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

proc = ob.ProcessFiles(
    dname=dname,
    options=options,
)

proc.get_darkcurrent()

proc.process_all()
data = proc.get_dataclass("Run00001").data

wavelengths = np.loadtxt('referenced_wavelengths.txt')

np.savez("data", delays=data['delays'], wavelengths=wavelengths, signal=data['mean'])


set_background = 2 # 0=exp, 1=mean, 2=linear
smooth_size = 3
butter_order = 4
t_echo = 42
t0_filter = 15 # ps
t0_filter_index = np.abs(data['delays'] - t0_filter).argmin()

plotlist = [1,12]
# 1. Full color delta plot

if 1 in plotlist:
    plt.figure(1)
    plt.pcolormesh(data['delays'], wavelengths, data['mean'].mean(2), cmap='RdBu')
    plt.clim(-0.01,0.01)
    plt.colorbar()
    for loc in np.arange(0, 500, t_echo):
        plt.axvline(loc, alpha=0.2, color='grey', linestyle=':', linewidth=1)
    plt.axhline(745, color='grey')
    plt.xlabel('delay [ps]')
    plt.ylabel('wavelength [nm]')
    plt.title('Reflection change')
    plt.show(block=False)

timestep = 0.1
delays = data['delays'][t0_filter_index:]
delta = data['mean'].mean(2)[:, t0_filter_index:]
delays_interp = np.arange(data['delays'][t0_filter_index], data['delays'][-1],timestep)
delta_interp = interpolate.interp1d(data['delays'], data['mean'].mean(2))(delays_interp) # starting after first peak

def exp_decay(x, a, b):
    return a * np.exp(-b * x)

if set_background == 0:
    bg_tag = 'Exp'
    param_list = []
    for i in range(2048):
        param, _ = curve_fit(exp_decay, delays_interp, delta_interp[i, :], bounds=([-0.01,0],[0.01,0.0008]))
        param_list.append(param)
        print(i, end='\r')
    delta_bgrm = np.array([delta[i, :] - exp_decay(delays, *param_list[i]) for i in range(2048)])
    delta_interp = np.array([delta_interp[i, :] - exp_decay(delays_interp, *param_list[i]) for i in range(2048)])

elif set_background == 1:
    bg_tag = 'Mean'
    delta_bgrm = delta - delta.mean(1)[:, np.newaxis]
    delta_interp = delta_interp - delta_interp.mean(1)[:, np.newaxis]

elif set_background == 2:
    bg_tag = 'Linear'
    mb_list = np.array([np.polyfit(delays_interp, delta_interp[i, :], 1) for i in range(2048)])
    delta_bgrm = np.array([delta[i, :] - (mb_list[i,0] * delays + mb_list[i,1]) for i in range(2048)])
    delta_interp = np.array([delta_interp[i, :] - (mb_list[i,0] * delays_interp + mb_list[i,1]) for i in range(2048)])

while True:
    hp_freq_GHz = input('HP cutoff frequency in GHz: ')
    hp_freq = float(hp_freq_GHz) * 1e9
    sos = signal.butter(butter_order, hp_freq, 'hp', fs=1e12, output='sos')

    single_wl = float(input('Trace single wavelength in nm: '))

    ## highpass filter
    delta_filtered = signal.sosfilt(sos, delta_bgrm, axis=1)

    # 2. Full color delta with HP
    if 2 in plotlist:
        plt.figure(2)
        plt.clf()
        plt.pcolormesh(delays, wavelengths, delta_filtered, cmap='RdBu')
        clim = max(np.abs([np.max(delta_filtered), np.min(delta_filtered)]))
        plt.clim(-clim,clim)
        plt.colorbar()
        plt.xlim(left=0)
        plt.xlabel('delay [ps]')
        plt.ylabel('wavelength [nm]')
        plt.title(f'Reflection change {hp_freq_GHz}GHz HP filter')
        plt.show(block=False)

    closest_wl_index = np.abs(wavelengths - single_wl).argmin()
    # delta trace of single wavelength
    delta_sgwl = delta[closest_wl_index, :]
    delta_sgwl_bgrm = delta_bgrm[closest_wl_index, :]
    delta_sgwl_filtered = delta_filtered[closest_wl_index, :]
    delta_sgwl_smooth = uniform_filter1d(delta_sgwl, size=smooth_size)
    delta_sgwl_bgrm_smooth = uniform_filter1d(delta_sgwl_bgrm, size=smooth_size)
    delta_sgwl_filtered_smooth = uniform_filter1d(delta_sgwl_filtered, size=smooth_size)
    
    # 3. delta traces at single wavelength
    if 3 in plotlist:
        plt.figure(3)
        plt.clf()
        plt.plot(delays, delta_sgwl_smooth, label='Original')
        plt.plot(delays, delta_sgwl_bgrm_smooth, label=f'{bg_tag} background removed')
        plt.plot(delays, delta_sgwl_filtered_smooth, label=f'{hp_freq_GHz} Ghz HP')
        plt.xlim(left=-20)
        plt.xlabel('delay [ps]')
        plt.ylabel(r'$\Delta R/R$')
        plt.title(f'Reflection change trace at {wavelengths[closest_wl_index]:.2f} nm')
        plt.legend()
        plt.show(block=False)

    outputlen = 2**11
    delta_filtered_fft = np.fft.rfft(delta_filtered, n=outputlen, axis=1)
    delta_filtered_fft_abs = np.abs(delta_filtered_fft)
    delta_filtered_fft_phase = np.angle(delta_filtered_fft)
    delta_fft = np.fft.rfft(delta_bgrm, n=outputlen, axis=1)
    delta_fft_abs = np.abs(delta_fft)
    delta_fft_phase = np.angle(delta_fft)
    delta_freqs = np.fft.rfftfreq(outputlen, d=1e-12) * 1e-9

    # 4. FFT with HP
    if 4 in plotlist:
        plt.figure(4)
        plt.clf()
        plt.plot(delta_freqs, delta_filtered_fft_abs[closest_wl_index, :])
        #plt.yscale('log')
        plt.xlim(0,40)
        plt.ylim(bottom=0)
        plt.xlabel('Frequency [GHz]')
        plt.ylabel('Amplitude [Arb. Unit]')
        plt.title(f'FFT of trace at {wavelengths[closest_wl_index]:.2f} nm {hp_freq_GHz}GHz HP')
        plt.show(block=False)

    # 5. FFT without HP
    if 5 in plotlist:
        plt.figure(5)
        plt.clf()
        plt.plot(delta_freqs, delta_fft_abs[closest_wl_index, :])
        #plt.yscale('log')
        plt.xlim(0,40)
        plt.ylim(bottom=0)
        plt.xlabel('Frequency [GHz]')
        plt.ylabel('Amplitude [Arb. Unit]')
        plt.title(f'FFT of trace at {wavelengths[closest_wl_index]:.2f} nm no HP')
        plt.show(block=False)
    
    # 6. 2D FFT with HP
    if 6 in plotlist:
        plt.figure(6)
        plt.clf()
        norm = 'linear'
        plt.pcolormesh(delta_freqs, wavelengths, delta_filtered_fft_abs, cmap='Blues', norm=norm)
        clim = np.max(delta_filtered_fft_abs[600:,:])
        plt.clim(0,clim)
        plt.colorbar()
        plt.xlim(0,40)
        plt.xlabel('frequency GHz')
        plt.ylabel('wavelengths [nm]')
        plt.title(f'FFT {hp_freq_GHz}GHz HP filter')
        plt.show(block=False)
    
    # 7. 2D FFT without HP
    if 7 in plotlist:
        plt.figure(7)
        plt.clf()
        norm = 'linear'
        plt.pcolormesh(delta_freqs, wavelengths, delta_fft_abs, cmap='Blues', norm=norm)
        clim = np.max(delta_fft_abs[600:,:])
        plt.clim(0,clim)
        plt.colorbar()
        plt.xlim(0,40)
        plt.xlabel('frequency GHz')
        plt.ylabel('wavelengths [nm]')
        plt.title(f'FFT no HP filter')
        plt.show(block=False)

    wl_list = np.array([single_wl - 30 + 10*i for i in range(7)])
    wl_indices = np.argmin(np.abs(wavelengths - wl_list[:, np.newaxis]), axis=1)

    # 8. FFT of several wavelengths
    if 8 in plotlist:
        plt.figure(8)
        plt.clf()
        color = iter(plt.cm.rainbow(np.linspace(0, 1, 7)))
        offset = 0
        offset_spacing = np.max(delta_fft_abs[wl_indices[4]])*1.1
        offset_list = np.arange(0, 7 * offset_spacing, offset_spacing)
        for i in wl_indices:
            plt.plot(delta_freqs, delta_fft_abs[i,:] + offset, c=next(color))
            offset += offset_spacing 
        plt.yticks(offset_list, color='None')
        plt.grid(axis='y')
        plt.xlim(0,40)
        plt.xlabel('Frequency (GHz)')
        plt.ylabel('Amplitude (Arb. Unit)')
        plt.title(f'FFT from multiple wavelengths without HP')
        plt.legend(wl_list, loc='lower right', reverse=True)
        plt.show(block=False)

    # 9. FFT of several wavelengths with HP
    if 9 in plotlist:
        plt.figure(9)
        plt.clf()
        color = iter(plt.cm.rainbow(np.linspace(0, 1, 7)))
        offset = 0
        offset_spacing = np.max(delta_filtered_fft_abs[wl_indices[4]])*1.1
        offset_list = np.arange(0, 7 * offset_spacing, offset_spacing)
        for i in wl_indices:
            plt.plot(delta_freqs, delta_filtered_fft_abs[i,:] + offset, c=next(color))
            offset += offset_spacing 
        plt.yticks(offset_list, color='None')
        plt.grid(axis='y')
        plt.xlim(0,40)
        plt.xlabel('Frequency (GHz)')
        plt.ylabel('Amplitude (Arb. Unit)')
        plt.title(f'FFT from multiple wavelengths with {hp_freq_GHz}GHz HP')
        plt.legend(wl_list, loc='lower right', reverse=True)
        plt.show(block=False)

    freq_list = np.array([2.5, 8, 25.5])
    freq_indices = np.argmin(np.abs(delta_freqs - freq_list[:, np.newaxis]), axis=1)
    
    # 10. Amplitudes at peak frequencies
    if 10 in plotlist:
        plt.figure(10)
        plt.clf()
        for i in freq_indices:
            plt.plot(wavelengths, uniform_filter1d(delta_fft_abs[:, i], size=12))
        plt.title('Amplitude of frequencies')
        plt.xlabel('wavelengths')
        plt.ylabel('Amplitudes Arb. Unit')
        plt.legend(freq_list)
        plt.show(block=False)

    # 11. Phase at peak frequencies
    if 11 in plotlist:
        plt.figure(11)
        plt.clf()
        plt.plot(wavelengths, uniform_filter1d(delta_fft_phase[:, freq_indices[0]], size=12)/np.pi)
        plt.plot(wavelengths, uniform_filter1d(delta_fft_phase[:, freq_indices[1]], size=12)/np.pi)
        plt.plot(wavelengths, uniform_filter1d(delta_fft_phase[:, freq_indices[2]], size=12)/np.pi)
        plt.ylim(-1,1)
        plt.title('Phase of frequencies')
        plt.xlabel('wavelengths')
        plt.ylabel(r'Phase [$\pi$]')
        plt.legend(freq_list)
        plt.show(block=False)

    wl_list = np.array([675,700,725,750,775,800,825])
    wl_indices = np.argmin(np.abs(wavelengths - wl_list[:, np.newaxis]), axis=1)

    # 12. Trace of several wavelengths without HP
    if 12 in plotlist:
        plt.figure(8)
        plt.clf()
        color = iter(plt.cm.rainbow(np.linspace(0, 1, 7)))
        offset = 0
        offset_spacing = np.max(delta_bgrm[wl_indices[4]])*1.1
        offset_list = np.arange(0, 7 * offset_spacing, offset_spacing)
        for i in wl_indices:
            plt.plot(delays, uniform_filter1d(delta_bgrm[i,:], size=3) + offset, c=next(color))
            offset += offset_spacing 
        plt.yticks(offset_list, color='None')
        plt.grid(axis='y')
        plt.xlabel('Delay (ps)')
        plt.ylabel(r'\Delta R/R_0')
        plt.title(f'Delta trace from multiple wavelengths without HP')
        plt.legend(wl_list, loc='lower right', reverse=True)
        plt.show(block=False)
    print()

print(".")

# %%
