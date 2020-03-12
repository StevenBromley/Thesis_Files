"""
Steve Bromley     02-13-19
Searching for peaks, frame by frame, in 6 cm data
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy import stats
import peakutils
import time
import csv
import glob
import re
from numpy import genfromtxt
start=time.time()

def read_prin(fil):
    f = open(fil)
    data = np.loadtxt(f,delimiter=',')
    inten = np.reshape(data[:,1],(-1,1340))
    return np.vstack([data[0:1340,0],inten])

def stdev(y):
   return(np.std(y))
   
def quad_std(x):
    sum = 0
    for i in range(0,len(x)):
        sum = sum + (x[i])**2
    mean = np.sqrt(sum)
    return(mean)   
#Max/min lambda for binning later; bounds will be adjusted depending on input
global_max_lambda = 0
global_min_lambda = 1000
#Save main working directory:
cwd = os.getcwd()
data_path = 'C:/Users/steve_000/Desktop/Research_Files/PRA_gold_spectroscopy/raw_data_noxrays_nosum/6cm'
os.chdir(data_path)   
#%%
#%%
manual_input = glob.glob("*.csv")
#Loop over ALL data files of interest in "manual_input" variable:
for data_in in manual_input:
    os.chdir(data_path)   
    pin =  genfromtxt(data_in,delimiter = ',')
    #Total Number of Frames in data file:
    frames = len(pin[:,0])
    #Final peak output array from gaussian fits
    all_peaks_out = np.empty((1340,200,2))
    #Select range of frames for analysis:
    if "181128" in data_in:
        frame1 = 6
        frame2 = 12
    if "181127" in data_in:
        frame1 = 62
        frame2 = 79
    #probe out frame range:
    #if "181126" in data_in:
     #   frame1 = 19
      #  frame2 = 23
   #Different trigger;for probe in, use below.
   #     frame1 = 4
    #    frame2 = 10
    if "181129" in data_in:
        frame1 = 55
        frame2 = 108
    if "181130" in data_in:
        frame1 = 45
        frame2 = 120
    #filter x-rays
    for frame in range(frame1,frame2):#range(frame1,frame2):
        baseline = peakutils.baseline(pin[frame,:], deg=3, max_it=100, tol=1e-5)
        indexes = peakutils.peak.indexes(pin[frame,:], thres=baseline + 100, min_dist=2, thres_abs='true')
        interp_peaks = peakutils.peak.interpolate(pin[0,:],pin[frame,:],indexes,width = 1)
        
    #Save output in 2 column format:
    #Peak Wavelength   Peak Counts
        numPeaks = len(indexes[:])
        all_peaks_out[0,frame,0] = frame
    #Get local max/min lambda from data and adjust globals:
        min_lambda = pin[0,0]
        max_lambda = pin[0,1339]   
        if min_lambda < global_min_lambda:
            global_min_lambda = min_lambda
        if max_lambda > global_max_lambda:
            global_max_lambda = max_lambda     
        for i in range(0,numPeaks):
            all_peaks_out[i+1,frame,0] = interp_peaks[i]
            all_peaks_out[i+1,frame,1] = pin[frame,indexes[i]] - baseline[i]
            if (interp_peaks[i] > max_lambda or interp_peaks[i] < min_lambda):
                    all_peaks_out[i+1,frame,0] = 0
    #For each row in output, save wavelength and counts in 2 column format         
    output = all_peaks_out[0:numPeaks,frame1:frame2,:]
    output = np.transpose(output)   
    #Array for storing 2d data used in next loop:
    output_2d = np.zeros((0,2))
    for j in range(0,len(output[0,:,0])):    
    #For a given row (j), add in all values from 0 -> end of frames data
        output_temp = np.zeros((len(output[0,0,1:]),2))  
        for k in range(0,len(output[0,0,1:])):
            output_temp[k,0] = output[0,j,k+1]
            output_temp[k,1] = output[1,j,k+1]
        output_2d = np.append(output_2d,output_temp, axis=0)
    #Sort:
    output_2d= output_2d[np.argsort(output_2d[:,0])]
    #Remove nan:
    output_2d = output_2d[output_2d.all(1)]
    #38.2 nm / 764 bins = .05nm bins; I used to use 390 
    num_Bins = 764
    #Mean peak positions:
    peaks_binned = scipy.stats.binned_statistic(output_2d[:,0],output_2d[:,0], statistic='mean', bins = num_Bins , range=(min_lambda,max_lambda))
    #Number of counts (intensity) in each bin
    counts_binned = scipy.stats.binned_statistic(output_2d[:,0],output_2d[:,1], statistic='sum', bins = num_Bins , range=(min_lambda,max_lambda))
    #Normalize all inputs to 1.25ms integration time (OH data already has this)
    if ("181126" or "181128" in data_in):
        counts_binned[0][:] = counts_binned[0][:] / 8
    #Number of peaks used to calculate the positions and counts in each bin
    pks_per_bin = scipy.stats.binned_statistic(output_2d[:,0],output_2d[:,1], statistic='count', bins = num_Bins , range=(min_lambda,max_lambda))
    #Calculate FINAL stdev (add stdev in quadrature):
    peaks_stdev = scipy.stats.binned_statistic(output_2d[:,0],output_2d[:,0], statistic=stdev, bins = num_Bins, range= (min_lambda,max_lambda))
    os.chdir(cwd)
    ##Export the avg peak and count to csv file
    a = peaks_binned[0][:]
    b = counts_binned[0][:]
    c = pks_per_bin[0][:]
    d = peaks_stdev[0][:]
    #Get shot # to save to output
    shot_num = re.findall("(\d+)",data_in)
    shot_num_int = int(shot_num[0])
    #Ex: output will be 18113026_500_540.csv -> 500- 540nm, shot num ...
    with open('{:3}_{:3}_{:3}.csv'.format(shot_num_int,min_lambda,max_lambda), 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(zip(a,b,c,d))
    f.close()
#%%
#%%
###Import all output CSV files and store in 2-D, sorted array
all_peak_data = np.empty((0,4))
data_files = glob.glob('1811*')
#Remove NaN entries from binning carried out previously and then sort
for data_file in data_files:    
    data_in = genfromtxt(data_file, delimiter = '\t')
    data_in = data_in[~np.isnan(data_in).any(axis=1)]
    all_peak_data = np.append(all_peak_data,data_in, axis=0)   
all_peak_data_sorted = all_peak_data[np.argsort(all_peak_data[:,0])]

###Rebinning the peaked data for all files
bin_size = .05
num_Bins = round((global_max_lambda - global_min_lambda)/bin_size)
testing2 = scipy.stats.binned_statistic(all_peak_data_sorted[:,0],all_peak_data_sorted[:,0],statistic='mean', bins = num_Bins, range= (global_min_lambda,global_max_lambda))
testing_count2 = scipy.stats.binned_statistic(all_peak_data_sorted[:,0],all_peak_data_sorted[:,1], statistic='sum', bins = num_Bins, range = (global_min_lambda, global_max_lambda))
testing2_pks_per_bin = scipy.stats.binned_statistic(all_peak_data_sorted[:,0],all_peak_data_sorted[:,2], statistic='sum', bins = num_Bins , range=(global_min_lambda,global_max_lambda))
testing_stdev = scipy.stats.binned_statistic(all_peak_data_sorted[:,0],all_peak_data_sorted[:,3],statistic=quad_std, bins = num_Bins, range= (global_min_lambda,global_max_lambda))
#%%
a1 = np.asarray(testing2[0][:])
b1 = np.asarray(testing_count2[0][:])
c1 = np.asarray(testing2_pks_per_bin[0][:])
d1 = np.asarray(testing_stdev[0][:])
abc1 = np.vstack((a1,b1,c1,d1))
#%%
final_lines_out = np.transpose(abc1)
final_lines_out = final_lines_out[~np.isnan(final_lines_out).any(axis=1)]
std_min = 0
final = np.empty((0,4))
final = np.vstack((final,final_lines_out[final_lines_out[:,3] > std_min]))
#Output FINAL line positions        
np.savetxt('final_line_list.txt',final,fmt='%.7f',delimiter=',')
#%%
#Remove temp files generated above
files_to_delete = glob.glob("1811*")
for f in files_to_delete:
    os.remove(f)