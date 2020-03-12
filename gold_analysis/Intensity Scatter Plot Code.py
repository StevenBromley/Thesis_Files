"""
Intensity Scatter Plot Script
SJB
Purpose: :ine intensities as a function of frame across all discharges
All defs and imports in steve.py
Required File Inputs:
-"shots.csv"; list of shot numbers
-A line list, either an import from a file OR a user-defined list of wavelengths, ex:
        peak_list_initial = np.array([400,500,600])
        peak_list_initial = genfromtxt('H lines.csv',delimiter =',')

If using on a PC, please make sure the datafile directories (below) are set properly for your PC!
"""

from steve import *
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

peak_list_initial = genfromtxt('Au_I_matches_07-15-19.txt',delimiter ='\t', usecols = 0, dtype = float)#np.array([486.14])
#gold list:
peak_pics = 'C:/Users/steve_000/Desktop/Research_Files/PRA_gold_spectroscopy/07-15-19/Au_integrated_intens_plots'
#%%
frame1 = 55
frame2 = 110
#load in the shots: the format here is from L - R go 0cm - 6cm (2nd probe)
#and going down increases the wavelength range
shot_list = genfromtxt('shots.csv',delimiter = ',')
#%%
for countr in range(0,17):
    #LOAD ALL OF THE PROBE DATA FOR THE WAVELENGTH RANGE
    #Load the data; 2nd probe:
    os.chdir(probe2)
    probe2_6cm = str(int(shot_list[countr,4])) + '.csv'
    probe2_6cm_dat = genfromtxt(probe2_6cm,delimiter =',')
    probe2_6cm_shot = str(int(shot_list[countr,4]))
    probe2_6cm_num = int(probe2_6cm_shot)
    
    min_lambda = probe2_6cm_dat[0,0]
    max_lambda = probe2_6cm_dat[0,1339]
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
        probe1_9cm_dat = genfromtxt(probe1_9cm,delimiter =',')
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
    ####Finished Data Load-in
    min_lambda = probe2_6cm_dat[0,0]
    max_lambda = probe2_6cm_dat[0,1339]
    peak_list = np.empty((0))
    #make the peak_list
    for i in range(0,len(peak_list_initial[:])):
        if (min_lambda < peak_list_initial[i] < max_lambda):
            peak_list = np.append(peak_list,peak_list_initial[i])
    #Save time by not peak finding if no peak is within the wavelength range:
    if (len(peak_list)==0):
        continue
    #%%Commence peak finding:    
    #This returns a 3d array of ALL peaks found between frames 1 and 2 in the data file; ex: next line
    #calculates all of the peaks between frame1 and frame2 for the shot data in probe1_0cm_dat
    #Repeat for all shots in a given wavelength range:
    pk_find_probe1_0cm = peak_find(frame1,frame2,probe1_0cm_dat)  
    pk_find_probe1_3cm = peak_find(frame1,frame2,probe1_3cm_dat)
    if (shot_list[countr,2] != 0 ):
        pk_find_probe1_6cm = peak_find(frame1,frame2,probe1_6cm_dat)    
        pk_find_probe1_9cm = peak_find(frame1,frame2,probe1_9cm_dat)    
    pk_find_probe2_6cm = peak_find(frame1,frame2,probe2_6cm_dat)    
    if (shot_list[countr,5] != 0 ):
        pk_find_probe_ni_0cm = peak_find(frame1,frame2,probe_ni_0cm_dat)    
    pk_find_probe_ni_3cm = peak_find(frame1,frame2,probe_ni_3cm_dat)    
#%%   
    for i in range(0,len(peak_list)):

    #generate the necessary data for scatter plot:
    #The scatter_data function returns a tuple where the data we want is in the component [1]
    #The function takes a peak input and searches through all of the data in the peak finding carried out above
    #and finds the closet peak. If the peak wavelength is within .07nm, i.e.~ +/- ~1 pixel 
    #then we consider it the same peak. It does this for the range of frames specified by frame1 to frame2
    #and then plots the peaks for all probe depths as a function of frame
        scatter_probe2_6cm = scatter_data(probe2_6cm_dat,pk_find_probe2_6cm[1],peak_list[i],frame1,frame2)
        scatter_probe1_0cm = scatter_data(probe1_0cm_dat,pk_find_probe1_0cm[1],peak_list[i],frame1,frame2)
        scatter_probe1_3cm = scatter_data(probe1_3cm_dat,pk_find_probe1_3cm[1],peak_list[i],frame1,frame2)
        if (shot_list[countr,2] != 0 ):
            scatter_probe1_6cm = scatter_data(probe1_6cm_dat,pk_find_probe1_6cm[1],peak_list[i],frame1,frame2)
            scatter_probe1_9cm = scatter_data(probe1_9cm_dat,pk_find_probe1_9cm[1],peak_list[i],frame1,frame2)
        if (shot_list[countr,5] != 0 ):
            scatter_probe_ni_0cm = scatter_data(probe_ni_0cm_dat,pk_find_probe_ni_0cm[1],peak_list[i],frame1,frame2)
        scatter_probe_ni_3cm = scatter_data(probe_ni_3cm_dat,pk_find_probe_ni_3cm[1],peak_list[i],frame1,frame2)
        os.chdir(peak_pics)
       #PLOT The Position Scatter Plots!
        plt.clf()
        plt.rcParams['font.size'] = 10
        plt.xlim(frame1,frame2)
        #Plot the peaks as a function of frame:
        plt.plot(scatter_probe1_0cm[:,3],scatter_probe1_0cm[:,2], linestyle='--', marker='v', color='k', label = 'Au Out')
        plt.plot(scatter_probe1_3cm[:,3],scatter_probe1_3cm[:,2], linestyle='--', marker='^', color='k', label = 'Au 3cm')
        if (shot_list[countr,2] != 0 ):
            plt.plot(scatter_probe1_6cm[:,3],scatter_probe1_6cm[:,2], linestyle='--', marker='s', color='k',label = 'Au 6cm')
            plt.plot(scatter_probe1_9cm[:,3],scatter_probe1_9cm[:,2], linestyle='--', marker='P', color='k', label = 'Au 9cm')
        #else:
            #plt.plot(0,0, label = 'Au(1) 6cm : N/A')
            #plt.plot(0,0, label = 'Au(1) 9cm: N/A')
        if (shot_list[countr,5] != 0 ):
            plt.plot(scatter_probe_ni_0cm[:,3],scatter_probe_ni_0cm[:,2], linestyle='--', marker='o', color='b',label = 'Ni Out')
        else:
            plt.plot(0,0,label = 'Ni: 0cm N/A')
        plt.plot(scatter_probe_ni_3cm[:,3],scatter_probe_ni_3cm[:,2], linestyle='--', marker='+', color='b',label = 'Ni 3cm')
        plt.plot(scatter_probe2_6cm[:,3],scatter_probe2_6cm[:,2], linestyle='--', marker='X', color='r',label = 'Au In 6m (Hot)')

            
        plt.legend(loc = 'upper right')
        plt.xlabel('Frame')
        plt.ylabel('Intensity (arb. units)')
        plt.title('peak_plot: {:} nm'.format(str(round(peak_list[i],2))))
        plt.savefig('peak_plot {:} nm.png'.format(peak_list[i],min_lambda,max_lambda))
os.chdir(cwd)

print (time.time()-start)