# %%
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

data = np.load('/Users/jeffreysuen/Photoacoustics Research/Data/grating whitelight measurement/Measurements/data.npz')
delays = data['delays']
wavelengths = data['wavelengths']
sig = data['signal'] * 100# (pixels, delays, 4 runs)
sig_mean = sig.mean(2) # in percentage

set_background = 2 # 0=exp, 1=mean, 2=linear
smooth_size = 3
butter_order = 4
t_echo = 130 * 2 / 6420 * 1e3
t0 = 2 # ps
t0_i = 58

deltacmap = 'seismic'

savefig = 1
plotlist = [6]
if 0 in plotlist:
    plt.figure(0)
    plt.plot(wavelengths, sig_mean[:,100])
# 1. Full color delta plot

if 1 in plotlist:
    plt.figure(1, figsize=(6.4,4.8))
    plt.pcolormesh(delays, wavelengths, sig_mean, cmap=deltacmap, rasterized=True)
    plt.clim(-1,1)
    plt.colorbar(label=r'$\Delta R/R$(%)')
    #for loc in np.arange(0, 500, t_echo):
    #    plt.axvline(loc, alpha=0.2, color='grey', linestyle=':', linewidth=2)
    plt.axhline(750, color='grey', linestyle='--', linewidth=0.5, alpha=0.5)
    plt.xlabel('delay (ps)')
    plt.ylabel('wavelength (nm)')
    plt.title('Reflection change')
    plt.tight_layout()
    plt.show(block=False)
    
    '''
    plt.figure(2, figsize=(6.4,4.8))
    plt.pcolormesh(delays, wavelengths, signal[:,:,0], cmap='RdBu', rasterized=True)
    plt.clim(-1,1)
    plt.colorbar(label=r'$\Delta R/R$(%)')
    #for loc in np.arange(0, 500, t_echo):
    #    plt.axvline(loc, alpha=0.2, color='grey', linestyle=':', linewidth=2)
    plt.axhline(750, color='grey', linestyle='--', linewidth=0.5, alpha=0.5)
    plt.xlabel('delay (ps)')
    plt.ylabel('wavelength (nm)')
    plt.title('Reflection change')
    plt.show(block=False)
    
    plt.figure(3, figsize=(6.4,4.8))
    plt.pcolormesh(delays, wavelengths, signal[:,:,1], cmap='RdBu', rasterized=True)
    plt.clim(-1,1)
    plt.colorbar(label=r'$\Delta R/R$(%)')
    #for loc in np.arange(0, 500, t_echo):
    #    plt.axvline(loc, alpha=0.2, color='grey', linestyle=':', linewidth=2)
    plt.axhline(750, color='grey', linestyle='--', linewidth=0.5, alpha=0.5)
    plt.xlabel('delay (ps)')
    plt.ylabel('wavelength (nm)')
    plt.title('Reflection change')
    plt.show(block=False)
    
    plt.figure(4, figsize=(6.4,4.8))
    plt.pcolormesh(delays, wavelengths, signal[:,:,2], cmap='RdBu', rasterized=True)
    plt.clim(-1,1)
    plt.colorbar(label=r'$\Delta R/R$(%)')
    #for loc in np.arange(0, 500, t_echo):
    #    plt.axvline(loc, alpha=0.2, color='grey', linestyle=':', linewidth=2)
    plt.axhline(750, color='grey', linestyle='--', linewidth=0.5, alpha=0.5)
    plt.xlabel('delay (ps)')
    plt.ylabel('wavelength (nm)')
    plt.title('Reflection change')
    plt.show(block=False)
    
    plt.figure(5, figsize=(6.4,4.8))
    plt.pcolormesh(delays, wavelengths, signal[:,:,3], cmap='RdBu', rasterized=True)
    plt.clim(-1,1)
    plt.colorbar(label=r'$\Delta R/R$(%)')
    #for loc in np.arange(0, 500, t_echo):
    #    plt.axvline(loc, alpha=0.2, color='grey', linestyle=':', linewidth=2)
    plt.axhline(750, color='grey', linestyle='--', linewidth=0.5, alpha=0.5)
    plt.xlabel('delay (ps)')
    plt.ylabel('wavelength (nm)')
    plt.title('Reflection change')
    plt.show(block=False)
    '''

delays = data['delays'][t0_i:]
sig = sig_mean[:, t0_i:]

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
    mb_list = np.array([np.polyfit(delays, sig[i,:], 1) for i in range(2048)])
    delta_bgrm = np.array([sig[i, :] - (mb_list[i,0] * delays + mb_list[i,1]) for i in range(2048)])

while True:
    hp_freq_GHz = 4
    hp_freq = float(hp_freq_GHz) * 1e9
    butter_order = 4
    sosff = signal.butter(butter_order, hp_freq, 'hp', fs=1e12, output='sos')
    sosf = signal.butter(butter_order*2, hp_freq, 'hp', fs=1e12, output='sos')

    single_wl = 730

    ## highpass filter
    starting_time = 15 # ps, first peak cut
    t1_i = 13
    delta_filtfilt = signal.sosfiltfilt(sosff, delta_bgrm, axis=1)
    delta_filt = signal.sosfilt(sosf, delta_bgrm, axis=1)

    # 2. Full color delta with HP
    if 2 in plotlist:
        plt.figure(2)
        plt.clf()
        plt.pcolormesh(delays[t1_i:], wavelengths, delta_filtfilt, cmap=deltacmap, rasterized=True)
        clim = max(np.abs([np.max(delta_filtfilt), np.min(delta_filtfilt)]))
        plt.clim(-clim,clim)
        plt.colorbar(label=r'$\Delta R/R$(%)')
        plt.xlim(left=0)
        plt.xlabel('delay [ps]')
        plt.ylabel('wavelength [nm]')
        plt.title(f'Reflection change {hp_freq_GHz}GHz HP filter')
        plt.tight_layout()
        plt.show(block=False)
    
    # 13. delta mean removed without hp
    if 2 in plotlist:
        plt.figure(13)
        plt.clf()
        plt.pcolormesh(delays, wavelengths, delta_bgrm, cmap=deltacmap, rasterized=True)
        clim = max(np.abs([np.max(delta_bgrm), np.min(delta_bgrm)]))
        plt.clim(-clim,clim)
        plt.colorbar(label=r'$\Delta R/R$(%)')
        plt.xlim(left=0)
        plt.xlabel('delay [ps]')
        plt.ylabel('wavelength [nm]')
        plt.title(f'Reflection change mean removed by row')
        plt.tight_layout()
        plt.show(block=False)
        
    closest_wl_index = np.abs(wavelengths - single_wl).argmin()
    # delta trace of single wavelength
    delta_sgwl = sig[closest_wl_index, :]
    delta_sgwl_bgrm = delta_bgrm[closest_wl_index, :]
    delta_sgwl_filtered = delta_filtfilt[closest_wl_index, :]
    delta_sgwl_smooth = uniform_filter1d(delta_sgwl, size=smooth_size)
    delta_sgwl_bgrm_smooth = uniform_filter1d(delta_sgwl_bgrm, size=smooth_size)
    delta_sgwl_filtered_smooth = uniform_filter1d(delta_sgwl_filtered, size=smooth_size)
    
    # 3. delta traces at single wavelength
    if 3 in plotlist:
        plt.figure(3)
        plt.clf()
        plt.plot(delays, delta_sgwl_smooth, label='Original')
        plt.plot(delays, delta_sgwl_bgrm_smooth, label=f'{bg_tag} background removed')
        plt.plot(delays[t1_i:], delta_sgwl_filtered_smooth, label=f'{hp_freq_GHz} Ghz HP')
        plt.xlim(left=-20)
        plt.xlabel('delay [ps]')
        plt.ylabel(r'$\Delta R/R$')
        plt.title(f'Reflection change trace at {wavelengths[closest_wl_index]:.2f} nm')
        plt.legend()
        plt.tight_layout()
        plt.show(block=False)

    outputlen = 2**11
    delta_filtered_fft = np.fft.rfft(delta_filtfilt, n=outputlen, axis=1)
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
        plt.tight_layout()
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
        plt.tight_layout()
        plt.show(block=False)
    
    # 6. 2D FFT with HP
    if 6 in plotlist:
        plt.figure(6)
        plt.clf()
        norm = 'linear'
        plt.pcolormesh(delta_freqs, wavelengths, delta_filtered_fft_abs, cmap='rainbow', norm='symlog', rasterized=True)
        #clim = np.max(delta_filtered_fft_abs[600:,:])
        #plt.clim(0,clim)
        plt.colorbar()
        plt.xlim(0,40)
        plt.xlabel('frequency GHz')
        plt.ylabel('wavelengths [nm]')
        plt.title(f'FFT {hp_freq_GHz}GHz HP filter')
        plt.tight_layout()
        plt.show(block=False)
    
    # 7. 2D FFT without HP
    if 7 in plotlist:
        plt.figure(7)
        plt.clf()
        norm = 'linear'
        plt.pcolormesh(delta_freqs, wavelengths, delta_fft_abs, cmap='rainbow', norm='symlog', rasterized=True)
        #clim = np.max(delta_fft_abs[600:,:])
        #plt.clim(0,clim)
        plt.colorbar()
        plt.xlim(0,40)
        plt.xlabel('frequency GHz')
        plt.ylabel('wavelengths [nm]')
        plt.title(f'FFT no HP filter')
        plt.tight_layout()
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
        plt.tight_layout()
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
        plt.tight_layout()
        plt.show(block=False)

    freq_list = np.array([24.7])
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
        plt.tight_layout()
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
        plt.tight_layout()
        plt.show(block=False)

    wl_interval = 10
    wl_list = np.array([750 - wl_interval*3 + i * wl_interval for i in range(7)])
    wl_indices = np.argmin(np.abs(wavelengths - wl_list[:, np.newaxis]), axis=1)

    # 12. Trace of several wavelengths without HP
    if 12 in plotlist:
        plt.figure(12, figsize=(3.1,4.8))
        color = iter(plt.cm.rainbow(np.linspace(0, 1, 7)))
        offset = 0
        offset_spacing = 0.4
        offset_list = np.arange(0, 7 * offset_spacing, offset_spacing)
        for i in wl_indices:
            plt.plot(delays, uniform_filter1d(delta_bgrm[i,:], size=5) + offset, c=next(color))
            offset += offset_spacing 
        plt.xticks([0,100,200,300,400,500])
        plt.yticks(offset_list)
        plt.grid(True, alpha=0.3)
        plt.ylim(-0.2,2.8)
        plt.xlabel('Delay (ps)')
        plt.ylabel(r'$\Delta R/R$ (%)')
        plt.title(f'mult wl no hp')
        plt.legend(wl_list, loc='lower right', reverse=True)
        plt.tight_layout()
        plt.show(block=False)
    
    if savefig==1:
        figdir = '/Users/jeffreysuen/Photoacoustics Research/Figures/Whitelight pump probe/pdf figs'
        for fignum in plt.get_fignums():
            if fignum != 1:
                plt.figure(fignum).set_size_inches(4,3.2)
                plt.tight_layout()  
            plt.savefig(f'{figdir}/figure {fignum}.pdf')
    print()
    
    

print(".")

# %%
