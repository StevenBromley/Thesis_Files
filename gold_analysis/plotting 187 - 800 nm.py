# -*- coding: utf-8 -*-
"""
SJB
Plotting 187 - 800 nm in one go
12-6-19
"""
import sys
sys.path.insert(0, 'C:/Users/steve_000/Desktop/Research_Files/PRA_gold_spectroscopy') 
from steve_lib import *
c_files = 'C:/Users/steve_000/Desktop/Research_Files/PRA_gold_spectroscopy/181130_params'
probe2 = 'C:/Users/steve_000/Desktop/Research_Files/PRA_gold_spectroscopy/raw_data_noxrays_nosum/2nd_probe'
regions = genfromtxt('regions_with_wavelengths.txt', delimiter ='\t')
shot_list = genfromtxt('shots.csv',delimiter = ',')
cwd = os.getcwd()
#%%
plt.clf()
plt.xlim(187,850)
plt.ylim(0,70000)

frame = 77
frame1 = frame
frame2 = frame

os.chdir(probe2)
for countr in range(0,17):
        probe2_6cm = str(int(shot_list[countr,4])) + '.csv'
        probe2_6cm_dat = genfromtxt(probe2_6cm,delimiter =',')
        probe2_6cm_shot = str(int(shot_list[countr,4]))
        probe2_6cm_num = int(probe2_6cm_shot)
        if (countr == 1):
            plt.plot(probe2_6cm_dat[0,:], probe2_6cm_dat[frame1,:])
        if (countr == 2):
            plt.plot(probe2_6cm_dat[0,:], probe2_6cm_dat[frame2,:])
        if ((countr != 1) and (countr != 2)):
            plt.plot(probe2_6cm_dat[0,:], probe2_6cm_dat[frame,:])
os.chdir(cwd)
    
plt.xlabel('Wavelength (nm)')
plt.ylabel('Intensity (arb. units)')


#Added on 1-06-20 to overlap Au lines over spectra
au0 = np.genfromtxt('au_I_lines_10-02-19.txt', delimiter ='\t', skip_header = 1)#, usecols =(0,1,2,3))
au1 = np.genfromtxt('au_II_lines_10-02-19.txt', delimiter ='\t', skip_header = 1)#, usecols =(0,1,2,3))
"""
for i in range(0,len(au0[:,0])):
    if (au0[i,0] > 0):
        plt.axvline(au0[i,1]/10, color = 'red')

for i in range(0,len(au1[:,0])):
    if (au1[i,0] > 0):
        plt.axvline(au1[i,1]/10, color = 'black', linestyle = '--')
"""