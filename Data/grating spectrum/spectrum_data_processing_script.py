# %%

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import oberon_package.dataprocessing as ob
import oberon_package.functions as fn
from oberon_package.parameters import Options
import os

abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

options = Options()

options.processing.ClearRM_measurement = "Run00005"
options.processing.AbsorptionRM_measurement = "Run00004"
options.processing.ref_shift = 2


proc = ob.ProcessFiles(dname, options)
Id = proc.get_darkcurrent()

r1 = ob.ReadPlus("Run00002", options, dname=dname) # grating
r1._read_txt()
r1.read_spectra_all()

r2 = ob.ReadPlus("Run00005", options, dname=dname) # flat
r2._read_txt()
r2.read_spectra_all()

I_grat = np.mean(r1.tmp, axis=0)[0,:]
I_flat = np.mean(r2.tmp, axis=0)[0,:]
Id_ref = Id[0,:]
reflection_spectrum = (I_grat - Id_ref)/(I_flat - Id_ref)

laser_wls = np.genfromtxt('referenced_wavelengths.txt')

Filenames = os.listdir(path='./18012024')
N_files = len(Filenames)

spectrometer_files = ['p_35_flat_Al_Subt2__0__15-06-27-274.txt',
                      'p_35_grat_Al_Subt2__0__15-07-34-975.txt',
                      'p_36_flat_Al_Subt2__0__15-10-02-275.txt', 
                      'p_36_grat_Al_Subt2__0__15-11-32-776.txt', 
                      'p_37_flat_Al_Subt2__0__15-13-43-877.txt', 
                      'p_37_grat_Al_Subt2__0__15-14-57-977.txt', 
                      'p_38_flat_Al_Subt2__0__15-17-42-878.txt', 
                      'p_38_grat_Al_Subt2__0__15-18-40-778.txt',] 
                    #   'p_39_flat_Al_Subt2__0__15-21-07-179.txt', 
                    #   'p_39_grat_Al_Subt2__0__15-21-59-479.txt']

raw_spectra_data = np.array([np.genfromtxt('./18012024/'+files, skip_header=14) for files in spectrometer_files])

gf_ratio = np.array([raw_spectra_data[i+1,534:911,1]/raw_spectra_data[i,534:911,1] for i in range(0,len(spectrometer_files),2)])
spectrometer_wls = raw_spectra_data[0,534:911,0]

plt.figure()
plt.plot(laser_wls, reflection_spectrum, c='black')
color = iter(cm.rainbow(np.linspace(0, 1, len(gf_ratio))))
for ratio in gf_ratio:
    plt.plot(spectrometer_wls, ratio, c=next(color))
plt.legend(['calibration', r'$35\degree$', r'$36\degree$', r'$37\degree$', r'$38\degree$'])
plt.xlabel('wavelength [nm]')
plt.ylabel('ratio')
plt.title('Reflection intensity ratio grating/flat Al')
plt.show(block=False)
    

print('.')