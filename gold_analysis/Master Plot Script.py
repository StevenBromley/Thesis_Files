"""
Position Scatter Plot Script
SJB
03-27-19
Purpose: Plot spectra for a single frame
All defs and imports in steve_lib.py
Required File Inputs:
-"shots.csv"; list of shot numbers
-A line list, either an import from a file OR a user-defined list of wavelengths, ex:
        peak_list_initial = np.array([400,500,600])
        peak_list_initial = genfromtxt('H lines.csv',delimiter =',')
"""
import sys
#define path for steve_lib.py
sys.path.insert(0, 'C:/Users/steve_000/Desktop/Research_Files/PRA_gold_spectroscopy')
#IMPORT EVERYTHING
from steve_lib import *
start=time.time()
#%%
cwd = os.getcwd()
#Directories where all of the raw shot files are stored:
probe2 = 'C:/Users/steve_000/Desktop/Research_Files/PRA_gold_spectroscopy/raw_data_noxrays_nosum/2nd_probe'
probe_out = 'C:/Users/steve_000/Desktop/Research_Files/PRA_gold_spectroscopy/raw_data_noxrays_nosum/probe_out'
probe3cm = 'C:/Users/steve_000/Desktop/Research_Files/PRA_gold_spectroscopy/raw_data_noxrays_nosum/3cm'
probe6cm = 'C:/Users/steve_000/Desktop/Research_Files/PRA_gold_spectroscopy/raw_data_noxrays_nosum/6cm'
probe9cm = 'C:/Users/steve_000/Desktop/Research_Files/PRA_gold_spectroscopy/raw_data_noxrays_nosum/9cm'
ni_0cm = 'C:/Users/steve_000/Desktop/Research_Files/PRA_gold_spectroscopy/nickel_raw_data_noxrays_nosum/ni_probe_out'
ni_3cm = 'C:/Users/steve_000/Desktop/Research_Files/PRA_gold_spectroscopy/nickel_raw_data_noxrays_nosum/ni_3cm'
ni_6cm = 'C:/Users/steve_000/Desktop/Research_Files/PRA_gold_spectroscopy/nickel_raw_data_noxrays_nosum/ni_6cm'
peak_list_initial = genfromtxt('fe1_matches_onlyreal3.txt',delimiter ='\t', usecols = 1)
shot_list = genfromtxt('shots.csv',delimiter = ',')
#%%
plt.clf()
countr = 6
frame = 81
frame_ni = 81
#Ex: Load the data; 2nd probe:
os.chdir(probe2)
probe2_6cm = str(int(shot_list[countr,4])) + '.csv'
probe2_6cm_dat = genfromtxt(probe2_6cm,delimiter =',')
probe2_6cm_shot = str(int(shot_list[countr,4]))
probe2_6cm_num = int(probe2_6cm_shot)
####Probe Out:
os.chdir(probe_out)
probe1_0cm = str(int(shot_list[countr,0])) + '.csv'
probe1_0cm_dat = genfromtxt(probe1_0cm,delimiter =',')
probe1_0cm_shot = str(int(shot_list[countr,0]))
probe1_0cm_num = int(probe1_0cm_shot)
###Probe1 3cm:
os.chdir(probe3cm)
probe1_3cm = str(int(shot_list[countr,1])) + '.csv'
probe1_3cm_dat = genfromtxt(probe1_3cm,delimiter =',')
probe1_3cm_shot = str(int(shot_list[countr,1]))
probe1_3cm_num = int(probe1_3cm_shot)
###Probe1 6cm:
os.chdir(probe6cm)
if (shot_list[countr,2] != 0 ):
    probe1_6cm = str(int(shot_list[countr,2])) + '.csv'
    probe1_6cm_dat = genfromtxt(probe1_6cm,delimiter =',')
    probe1_6cm_shot = str(int(shot_list[countr,2])) 
    probe1_6cm_num = int(probe1_6cm_shot)
    ###Probe1 9cm:
    os.chdir(probe9cm)
    probe1_9cm = str(int(shot_list[countr,3])) + '.csv'
    probe1_9cm_dat= genfromtxt(probe1_9cm,delimiter =',')
    probe1_9cm_shot = str(int(shot_list[countr,3]))
    probe1_9cm_num = int(probe1_9cm_shot)
#Nickel:    
if (shot_list[countr,5] != 0 ):
    os.chdir(ni_0cm)
    probe_ni_0cm = str(int(shot_list[countr,5])) + '.csv'
    probe_ni_0cm_dat = genfromtxt(probe_ni_0cm,delimiter =',')
    probe_ni_0cm_shot = str(int(shot_list[countr,5]))
    probe_ni_0cm_num = int(probe_ni_0cm_shot)

os.chdir(ni_3cm)
probe_ni_3cm = str(int(shot_list[countr,6])) + '.csv'
probe_ni_3cm_dat = genfromtxt(probe_ni_3cm,delimiter =',')
probe_ni_3cm_shot = str(int(shot_list[countr,6]))
probe_ni_3cm_num = int(probe_ni_3cm_shot)
#%
#only plot if file exists, i.e. shot_list var not 0
if (shot_list[countr,7] != 0 ):
    os.chdir(ni_6cm)
    probe_ni_6cm = str(int(shot_list[countr,7])) + '.csv'
    probe_ni_6cm_dat = genfromtxt(probe_ni_6cm,delimiter =',')
    probe_ni_6cm_shot = str(int(shot_list[countr,7]))
    probe_ni_6cm_num = int(probe_ni_6cm_shot)

min_lambda = probe2_6cm_dat[0,0]
max_lambda = probe2_6cm_dat[0,1339]
data_in = probe2_6cm_dat
os.chdir(cwd)

plt.grid(True)
plt.plot(probe2_6cm_dat[0,:],probe2_6cm_dat[frame,:], color = 'r', label = 'Au$_2$ 6cm')
#plt.plot(probe1_0cm_dat[0,:],probe1_0cm_dat[frame,:], color = 'k', label = 'Au(1) 0')
#plt.plot(probe1_3cm_dat[0,:],probe1_3cm_dat[frame,:], color = 'b', label = 'Au(1) 3')
#plt.plot(probe1_6cm_dat[0,:],probe1_6cm_dat[frame,:], color = 'g', label = 'Au(1) 6')
#plt.plot(probe_ni_0cm_dat[0,:],probe_ni_0cm_dat[frame,:], color = 'orange', label = 'Ni(1) 0')
plt.plot(probe_ni_3cm_dat[0,:],probe_ni_3cm_dat[frame_ni,:], color = 'purple', label = 'Ni 3cm')

if (shot_list[countr,7] != 0 ):
    plt.plot(probe_ni_6cm_dat[0,:],probe_ni_6cm_dat[frame_ni,:], color = 'orange', label = 'Ni 6cm')

plt.xlim(min_lambda,max_lambda)
plt.ylim(0,70000)
plt.xlabel('Wavelength (nm)')
plt.ylabel('Intensity (arb. units)')
plt.legend()
for i in range(0,len(peak_list_initial)):
    plt.axvline(peak_list_initial[i],linestyle = '--', color = 'k')

#for i in range(0,len(peak_list_initial2)):
#    plt.axvline(peak_list_initial2[i],linestyle = '--', color = 'r')
            
#for i in range(0,len(fe0_LTE[:,0])):
#    plt.axvline(fe0_LTE[i,0],fe0_LTE[0,1]*scale, linewidth = .7, ymax = fe0_LTE[0,1]*scale,label = 'FeI LTE',color ='blue')

plt.plot()
print (time.time()-start)