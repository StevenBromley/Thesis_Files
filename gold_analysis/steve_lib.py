"""
steve.py
SJB 
This is a file containing all of the imports needed to run my analysis codes for the CTH data.
This is the master file of my libraries. To import to a script, use the following code:
    
import sys
#define path for steve_lib.py
sys.path.insert(0, 'C:/Users/steve_000/Desktop/Research_Files/PRA_gold_spectroscopy') 
from steve_lib import *



"""
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy import stats
from scipy.optimize import curve_fit
from numpy import inf
import peakutils
import time
import csv
import glob
import re
from numpy import genfromtxt
start=time.time()
from scipy.integrate import *
from math import erf
from math import *
######
"""Baseline estimation algorithms."""

import numpy as np
import scipy.linalg as LA
import math

#####################
start=time.time()
#%%
def TTF(val, start):
    #function to convert a CTH time to a frame
    #start: spectrometer trigger start time in seconds
    #val: time (in second) to convert
    val2 = 1000 * (val - start) / 1.25
    return val2

#%%
#Intensity function version 3 (01-06-19)
#Incorporates calibration from peak overlap regions

def intensity_search_v3(array,frame_in,normalize,scale_max,method,regions,data_directory,frame_source,calib,correct):
    ###########################################################################################################################################
    ######################  Explanation of the parameters   ###################################################################################
    #                                                                                                                                         # 
    #   array: input array of line positions                                                                                                  # 
    #   frame: frame you want the intensities from                                                                                            # 
    #   normalize: boolean; if 'True', normalize the largest value to scale_max                                                               #
    #   scale_max: absolute maximum by which the found line intensities are normalized                                                        #
    #   method: select either 'relative' or 'absolute'                                                                                        # 
    #         'relative': normalizes to scale_max within each wavelength region                                                               #
    #         'absolute': normalizes to scale_max over all wavelength ranges (not reliable)                                                   #
    #   regions: array containing the wavelength regions from the raw data. Ex:                                                               #
    #            Format: (region #) (lower lambda) (upper lambda) (frame for intensities) (deg. for baseline polynomial)                      #
    #   frame_source: 'user' uses the value of 'frame' in the function call, or 'file' uses the col #4 from regions array                     #
    #   calib: intensity response of detector on 0 - 100 scale; 2 col format: Wavelength (nm) | Response. (defined 0 as first element)        #
    #   correct: boolean controlling whether or not to re-calibration with the above calib file                                               #
    ###########################################################################################################################################
    cwd = os.getcwd()
    #define the data directory (probe2), output dir. for plots (frame plots), and shot_list (wav. regions)
    shot_list = genfromtxt('shots.csv',delimiter = ',')
    lines = array
    os.chdir(cwd)
    temp = np.empty((1))
    #initial output of intensity search is  intensity is intens_arr
    intens_arr = np.empty((0,4))
    for countr in range(0,17):
        #load the data for the shot in wavelength range defined by countr:
        os.chdir(data_directory)
        probe2_6cm = str(int(shot_list[countr,4])) + '.csv'
        probe2_6cm_dat = genfromtxt(probe2_6cm,delimiter =',')
        probe2_6cm_shot = str(int(shot_list[countr,4]))
        probe2_6cm_num = int(probe2_6cm_shot)
        data_in = probe2_6cm_dat
        x_min = data_in[0,0]
        x_max = data_in[0,1339]
        lines_in_region = np.empty((0))
        #temp arr for holding lines
        temp1 = np.empty((1))
        #grab the lines within this wavelength window
        for i in range(0,len(lines)):    
            if (x_min < lines[i] < x_max):
                temp1 = lines[i]
                lines_in_region = np.append(lines_in_region,temp1)
        lines_in_region = lines_in_region[lines_in_region[:].argsort()]                
        #%Find the peaks using my standard peak finding function:
        #frame1 = 60
        #frame2 = 120
        #peak list from peak finding code. Do for all frames from frame1 -frame2
        if (frame_source == 'user'):
            deg = int(6)
            frame = int(frame_in)
        if (frame_source =='file'):
            deg = int(regions[countr,4])
            frame = int(regions[countr,3])
        peak_ls = peak_find(frame,frame+1,data_in)
        #loop over all frames :)
        #Find the baseline values & initialize values for the total gaussian and the spectra function
        total_gauss = 0
        #modified the peakutils baseline function to return the coefficients too :)
        baseline_tot = baseline_mod(data_in[frame,:], deg, max_it=100, tol=1e-5)
        baseline = baseline_tot[0]
        #define an array of coefficients for poly1d function:
        #fit our x values to a 'deg' degree polynomial for the baseline y values; return coeffs
        poly_params = np.polyfit(data_in[0,:], baseline, deg)
        #define a polynomial with the resulting coeffs. Use the function 'poly' below for integrating
        poly = np.poly1d(poly_params)
        poly_dat = poly(data_in[0,:])
        #the baseline acts as the starting point for constructing the 'synthetic spectra'
        synthetic_spec = baseline
        for i in range(0,len(lines_in_region)):
            #open up an output file for troubleshooting if need be:
            #file = open('testfile_frame_{:}.txt'.format(lines_in_region), 'w')
            #tempa array to hold each peaks data before saving:
            temp = np.empty((1,4))
            #only look for the peak & fit IF the peak is in this wavelength region:
            location = find_nearest_index(data_in[0,:], lines_in_region[i])
            for j in range(0,len(data_in[0,:])):
                #check to make sure the peak finding found the peak in this frame:
                if (abs(lines_in_region[i] - peak_ls[1][frame,j,0]) < .04):
                    #make sure the peak isn't at the borders:
                    if (5 < location < 1334):
                        #try to fit; otherwise throw an error and go to next peak:
                        try:
                            #Set the initial parameters for height, center, and width
                           #%%
                            #location = 1057
                            A_0 = data_in[frame,location]
                            #mu_0 = data_in[0,location]
                            mu_0 = peak_ls[1][frame,j,0]
                            sigma_0 = 0.1
                            p0 = [A_0,mu_0,sigma_0] #C_0]
                            #changed from .03 to .04 on 10-02-19
                            c1 = .04#c1 control the constraint on the line center
                            #limit the peak height unless the peak is saturated
                            if (data_in[frame,location] == 65535):
                                count_wiggle = 200000
                            else:
                                count_wiggle = 100                           
                            bounds_init = ((data_in[frame,location] - 100 , peak_ls[1][frame,j,0] - c1,0), \
                                           (data_in[frame,location] + count_wiggle, peak_ls[1][frame,j,0] + c1, 0.14))
                            #fit to a gaussian using the curve_fit module:
                            coeff, var_matrix = curve_fit(gaussian, data_in[0,:], data_in[frame,:], p0=p0, \
                                bounds = bounds_init) 
                            fit = gaussian(data_in[0,:],*coeff)
                            #plt.plot(data_in[0,:], fit)
                            #%%
                            #intersections between gaussian and baseline                    
                            intersections = np.argwhere(np.diff(np.sign(fit - baseline))).flatten()
                            #new x limits are the intersections between baseline and gaussian functions
                            lower_x = data_in[0,intersections[0]]
                            upper_x = data_in[0,intersections[1]]
                            #background counts for this peak:
                            bkg_counts = quad(poly,lower_x,upper_x)
                            #numerically integrated intens:
                            int_intens = quad(gaussian, lower_x , upper_x, args = (coeff[0],coeff[1],coeff[2]))
                            #analytically integrated intensity; save both later
                            analytic_intens = analytic_counts(lower_x, upper_x, coeff[0],coeff[1],coeff[2])
                            #%linspace integration; outdated but keeping for future if needed]
                            #x_data = data_in[0,:]
                            #x_range = np.linspace(data_in[0,location] - 0.5, data_in[0,location] + 0.5 ,num=10000, dtype = float)
                            #linspace_intens = quad(gaussian, x_range[0] , max(x_range), args = (coeff[0],coeff[1],coeff[2]))
                            #print('Numerical Intens: {:}'.format(int_intens[0] - bkg_counts[0]))
                            #print('Analytic Intens: {:}'.format(analytic_intens - bkg_counts[0]))
                            #print(coeff[1])
                            #print(coeff[2])
                            #%find the total gaussian fit function & calc. the synth. spectra
                            total_gauss = gauss_sum(total_gauss, data_in[0,:], coeff[0],coeff[1],coeff[2])
                            synthetic_spec = synthetic_spectra(synthetic_spec, fit[:],intersections)
                            #save wavelength, numerical intens, analytic intens, and frame # for each peak
                            temp[0,0] = lines_in_region[i]
                            #%%
                            temp[0,1] = int_intens[0] - bkg_counts[0]
                            if (temp[0,1] < 0):
                                temp[0,1] = 0
                            temp[0,2] = analytic_intens - bkg_counts[0]
                            temp[0,3] = frame 
                            intens_arr = np.append(intens_arr,temp,axis=0)
                        except RuntimeError:
                            print('Error - curve_fit failed for peak at {:}'.format(lines_in_region[i]))
                            print('Peak {:} added to output with intensity 0'.format(lines_in_region[i]))
                            
                            #Add in a 0 intensity entry for the failed peak
                            #Added on 1-7-19
                            temp = np.empty((1,4))
                            temp[0,0] = lines_in_region[i]
                            temp[0,1] = 0
                            temp[0,2] = 0
                            temp[0,3] = frame 
                            intens_arr = np.append(intens_arr,temp,axis=0)
                        except IndexError:
                            print('Error - No Intersections (peak {:} DNE)'.format(round(lines_in_region[i],5)))   
        #Commented out on 10-02-19 as a precautionary measure. Will re-enable if necessary
        """
            if (lines_in_region not in temp):
                #try to fit a peak in case the peak was not detected. If it works, save it, otherwise trash it
                try:
                    #initial guesses for the fit parameters:
                    location = find_nearest_index(data_in[0,:], lines_in_region[i])
                    A_0 = data_in[frame,location]
                    mu_0 = data_in[0,location]
                    sigma_0 = 0.1
                    p0 = [A_0,mu_0,sigma_0] #C_0]
                    c1 = .02#c1 control the constraint on the line center
                    count_wiggle = 100 
                    #wiggle room in the counts = count_wiggle. Creative variable name
                    bounds_init = ((data_in[frame,location] - 100 , data_in[0,location] - c1,0), \
                                           (data_in[frame,location] + count_wiggle, data_in[0,location] + c1, 0.14))                          
                    coeff, var_matrix = curve_fit(gaussian, data_in[0,:], data_in[frame,:], p0=p0, \
                                bounds = bounds_init) 
                    fit = gaussian(data_in[0,:],*coeff)       
                    #find where the baseline function & gaussian intersect:
                    intersections = np.argwhere(np.diff(np.sign(fit - baseline))).flatten()
                    #save the upper and lower bounds of the region containing the peak above the basline
                    lower_x = data_in[0,intersections[0]]
                    upper_x = data_in[0,intersections[1]]
                    #numerically integrate the background profile:
                    bkg_counts = quad(poly,lower_x,upper_x)
                    int_intens = quad(gaussian, lower_x , upper_x, args = (coeff[0],coeff[1],coeff[2]))
                    analytic_intens = analytic_counts(lower_x, upper_x, coeff[0],coeff[1],coeff[2])
                    total_gauss = gauss_sum(total_gauss, data_in[0,:], coeff[0],coeff[1],coeff[2])
                    synthetic_spec = synthetic_spectra(synthetic_spec, fit[:],intersections)
                    #save wavelength, numerical intens, analytic intens, and frame # for each peak
                    temp[0,0] = lines_in_region[i]
                    #save the integrated intensity minus the background contribution:
                    temp[0,1] = int_intens[0] - bkg_counts[0]
                    temp[0,2] = analytic_intens - bkg_counts[0]
                    temp[0,3] = frame 
                    intens_arr = np.append(intens_arr,temp,axis=0)
                except RuntimeError:
                    print('Error - curve_fit failed for peak at {:}'.format(lines_in_region[i]))
                except IndexError:
                    print('Error - No Intersections (peak {:} DNE)'.format(round(lines_in_region[i],5)))
        """
    intensity_out = np.empty((0,2))
    os.chdir(cwd)
    #sort the output and return the integrated intensities
    for i in range(0,len(lines[:])):
        temp_intens = np.empty((1,2))
        temp_intens[0,0] = lines[i]
        for j in range(0,len(intens_arr[:,0])):
            if (intens_arr[j,0] == lines[i]):
                temp_intens[0,1] = temp_intens[0,1] + intens_arr[j,1]
        intensity_out = np.append(intensity_out, temp_intens, axis=0)
   #%%
    #Added on 01-06-20 to incorporate spectrometer response calibration
    if ('correct' == True):
        for p in range(0,len(intensity_out[:,0])):
            #find the correction factor needed to apply the calibration
            #search in the calib (column 0) for nearest index of the entry which has the same
            #wavelength as the intensity_out array value
            calib_val = calib[find_nearest_index(calib[:,0],intensity_out[p,0]),1]
            intensity_out[p,1] = intensity_out[p,1] * 100 / calib_val
    
    #Is that it? That easy? Then just apply the absolute normalization below :)
    #%%
    #SJB 09-06-19 I know the 'True' and 'absolute' functions work.
    if ((normalize == 'True') and (method == 'absolute')):
        max_val = np.amax(intensity_out[:,1])
        intensity_out[:,1] = scale_max * intensity_out[:,1] / max_val  
    #work with a copy of intensity_out. We'll use this copy to change the values in intensity_out to their normalized values
    #later
    intensity_out_temp = intensity_out
    #Not working
    if ((normalize == 'True') and (method == 'relative')):
        for r in range(0,len(regions[:,0])):
            lower_x = regions[r,1]
            upper_x = regions[r,2]
            region_list = np.empty((0,2))
            #find the lines within each wavelength region:
            for i in range(0,len(intensity_out_temp[:,0])):
                temp = np.empty((1,2))
                if (lower_x < intensity_out_temp[i,0] < upper_x):
                    temp[0,0] = intensity_out_temp[i,0]
                    temp[0,1] = intensity_out_temp[i,1]
                    region_list = np.append(region_list,temp, axis = 0)
            
            if (r == 10):
                np.savetxt('region_list_test.txt', region_list)
            #find the maximum for lines in this wavelength region 
            try:
                max_val = np.amax(region_list[:,1])
                #normalize to the maximum value in this region
                region_list[:,1] = scale_max * region_list[:,1] / max_val
               
                #replace the intensities in the master list with the normalized values:
                for i in range(0,len(intensity_out_temp[:,0])):
                        j_max = len(region_list)    
                        for j in range(0,len(region_list)):
                            #Added on 10-15-19 to fix issue with blend intensity                                
                            if (region_list[j,0] == intensity_out_temp[i,0]):
                                if (len(region_list) < 2):
                                    intensity_out_temp[i,1] = region_list[j,1]
                                if (region_list[j,0] != region_list[j-1,0]):
                                    intensity_out_temp[i,1] = region_list[j,1]

            except ValueError:
                #STEVE I should add in check that if the peak IS NOT found then it adds a 0 to the array
                #will do this next time so the entire intensity code is self-contained
                #and won't rely on any post processing or checking
                print('No Peaks in wavelength range {:} : {:} - {:} nm'.format(countr,regions[r,1],regions[r,2]))
    intensity_out = intensity_out_temp
    return intensity_out





#%%
def one_peak_intens(line,data,peak_found_list,frame1,frame2):
    #######################################################################################################
    ######################  Explanation of the parameters   ###############################################
    #                                                                                                     #  
    #   line: a single line that you want the integrated intensity for                                    #
    #   data: this is the raw data array you will fit to; ex: a probe 2 shot file                         #
    #   peak_found_list: the output array of the the [1] component of the tuple returned                  #
    #                    from the peak_find function                                                      #              
    #   frame1: start of frame range over which to search                                                 # 
    #   frame2: end of frame range over which to search                                                   #  
    #   normalize: boolean; if 'True', normalize the largest value to scale_max                           #
    #   scale_max: absolute maximum by which the found line intensities are normalized                    #
    #                                                                                                     #
    #######################################################################################################

    cwd = os.getcwd()
    peak_ls = peak_found_list
    output = np.empty((0,4))
    for frame in range(frame1, frame2):
        deg = 10
        baseline_tot = baseline_mod(data[frame,:], deg, max_it=100, tol=1e-5)
        baseline = baseline_tot[0]
        #finds coefficients:
        poly_params = np.polyfit(data[0,:], baseline, deg)
        #define a polynomial with the resulting coeffs. Use the function 'poly' below for integrating
        poly = np.poly1d(poly_params)
        poly_dat = poly(data[0,:])
        location = find_nearest_index(data[0,:], line)
       
        try:
        #Set the initial parameters for height, center, and width
            A_0 = data[frame,location]
            mu_0 = data[0,location]
            sigma_0 = 0.1
            p0 = [A_0,mu_0,sigma_0] #C_0]
            c1 = .03#c1 control the constraint on the line center
            #limit the peak height unless the peak is saturated
            if (data[frame,location] == 65535):
                count_wiggle = 200000
            else:
                count_wiggle = 100                           
            bounds_init = ((data[frame,location] - 100 , data[0,location] - c1,0), \
                           (data[frame,location] + count_wiggle, data[0,location] + c1, 0.14))
            #fit to a gaussian using the curve_fit module:
            coeff, var_matrix = curve_fit(gaussian, data[0,:], data[frame,:], p0=p0, \
            bounds = bounds_init) 
            fit = gaussian(data[0,:],*coeff)
            #intersections between gaussian and baseline                    
            intersections = np.argwhere(np.diff(np.sign(fit - baseline))).flatten()
            #new x limits are the intersections between baseline and gaussian functions
            lower_x = data[0,intersections[0]]
            upper_x = data[0,intersections[1]]
            #background counts for this peak:
            bkg_counts = quad(poly,lower_x,upper_x)
            #numerically integrated intens:
            int_intens = quad(gaussian, lower_x , upper_x, args = (coeff[0],coeff[1],coeff[2]))
            #analytically integrated intensity; save both later
            analytic_intens = analytic_counts(lower_x, upper_x, coeff[0],coeff[1],coeff[2])
            #%linspace integration; outdated but keeping for future if needed]
            #x_data = data_in[0,:]
            #x_range = np.linspace(data_in[0,location] - 0.5, data_in[0,location] + 0.5 ,num=10000, dtype = float)
            #linspace_intens = quad(gaussian, x_range[0] , max(x_range), args = (coeff[0],coeff[1],coeff[2]))
            #print('Numerical Intens: {:}'.format(int_intens[0] - bkg_counts[0]))
            #print('Analytic Intens: {:}'.format(analytic_intens - bkg_counts[0]))
            #print(coeff[1])
            #print(coeff[2])
            #%find the total gaussian fit function & calc. the synth. spectra
            #save wavelength, numerical intens, analytic intens, and frame # for each peak
            output_temp = np.empty((1,4))
            output_temp[0,0] = line
            output_temp[0,1] = analytic_intens
            output_temp[0,2] = int_intens[0]
            output_temp[0,3] = frame
            output = np.append(output,output_temp,axis=0)
        except RuntimeError:
            #print('Error - curve_fit failed for peak at {:}'.format(line))
            #file.close()
            output = np.empty((1,4))
            output_temp[0,0] = line
            output_temp[0,1] = 0
            output_temp[0,2] = 0
            output_temp[0,3] = frame
            output = np.append(output,output_temp,axis=0)

        except IndexError:
            #print('Error - No Intersections (peak {:} DNE)'.format(round(line,5))) 
            output_temp = np.empty((1,4))
            output_temp[0,0] = line
            output_temp[0,1] = 0
            output_temp[0,2] = 0
            output_temp[0,3] = frame
            output = np.append(output,output_temp,axis=0)
                            
    return(output)

def intensity_search_v2(array,frame_in,normalize,scale_max,method,regions,data_directory,frame_source):
    ###########################################################################################################################################
    ######################  Explanation of the parameters   ###################################################################################
    #                                                                                                                                         # 
    #   array: input array of line positions                                                                                                  # 
    #   frame: frame you want the intensities from                                                                                            # 
    #   normalize: boolean; if 'True', normalize the largest value to scale_max                                                               #
    #   scale_max: absolute maximum by which the found line intensities are normalized                                                        #
    #   method: select either 'relative' or 'absolute'                                                                                        # 
    #         'relative': normalizes to scale_max within each wavelength region                                                               #
    #         'absolute': normalizes to scale_max over all wavelength ranges (not reliable)                                                   #
    #   regions: array containing the wavelength regions from the raw data. Ex:                                                               #
    #            Format: (region #) (lower lambda) (upper lambda) (frame for intensities) (deg. for baseline polynomial)                      #
    #   frame_source: 'user' uses the value of 'frame' in the function call, or 'file' uses the col #4 from regions array                     #
    ###########################################################################################################################################
    cwd = os.getcwd()
    #define the data directory (probe2), output dir. for plots (frame plots), and shot_list (wav. regions)
    shot_list = genfromtxt('shots.csv',delimiter = ',')
    lines = array
    os.chdir(cwd)
    temp = np.empty((1))
    #initial output of intensity search is  intensity is intens_arr
    intens_arr = np.empty((0,4))
    for countr in range(0,17):
        #load the data for the shot in wavelength range defined by countr:
        os.chdir(data_directory)
        probe2_6cm = str(int(shot_list[countr,4])) + '.csv'
        probe2_6cm_dat = genfromtxt(probe2_6cm,delimiter =',')
        probe2_6cm_shot = str(int(shot_list[countr,4]))
        probe2_6cm_num = int(probe2_6cm_shot)
        data_in = probe2_6cm_dat
        x_min = data_in[0,0]
        x_max = data_in[0,1339]
        lines_in_region = np.empty((0))
        #temp arr for holding lines
        temp1 = np.empty((1))
        #grab the lines within this wavelength window
        for i in range(0,len(lines)):    
            if (x_min < lines[i] < x_max):
                temp1 = lines[i]
                lines_in_region = np.append(lines_in_region,temp1)
        lines_in_region = lines_in_region[lines_in_region[:].argsort()]                
        #%Find the peaks using my standard peak finding function:
        #frame1 = 60
        #frame2 = 120
        #peak list from peak finding code. Do for all frames from frame1 -frame2
    
        if (frame_source == 'user'):
            deg = int(6)
            frame = int(frame_in)

        if (frame_source =='file'):
            deg = int(regions[countr,4])
            frame = int(regions[countr,3])
            
        
        
        peak_ls = peak_find(frame,frame+1,data_in)
        #loop over all frames :)
        #Find the baseline values & initialize values for the total gaussian and the spectra function
        total_gauss = 0
        #modified the peakutils baseline function to return the coefficients too :)
    
        baseline_tot = baseline_mod(data_in[frame,:], deg, max_it=100, tol=1e-5)
        baseline = baseline_tot[0]
        #define an array of coefficients for poly1d function:
        #fit our x values to a 'deg' degree polynomial for the baseline y values; return coeffs
        poly_params = np.polyfit(data_in[0,:], baseline, deg)
        #define a polynomial with the resulting coeffs. Use the function 'poly' below for integrating
        poly = np.poly1d(poly_params)
        poly_dat = poly(data_in[0,:])
        #the baseline acts as the starting point for constructing the 'synthetic spectra'
        synthetic_spec = baseline
        for i in range(0,len(lines_in_region)):
            #open up an output file for troubleshooting if need be:
            #file = open('testfile_frame_{:}.txt'.format(lines_in_region), 'w')
            #tempa array to hold each peaks data before saving:
            temp = np.empty((1,4))
            #only look for the peak & fit IF the peak is in this wavelength region:
            location = find_nearest_index(data_in[0,:], lines_in_region[i])
            for j in range(0,len(data_in[0,:])):
                #check to make sure the peak finding found the peak in this frame:
                if (abs(lines_in_region[i] - peak_ls[1][frame,j,0]) < .04):
                    #make sure the peak isn't at the borders:
                    if (5 < location < 1334):
                        #try to fit; otherwise throw an error and go to next peak:
                        try:
                            #Set the initial parameters for height, center, and width
                           #%%
                            #location = 1057
                            A_0 = data_in[frame,location]
                            #mu_0 = data_in[0,location]
                            mu_0 = peak_ls[1][frame,j,0]
                            sigma_0 = 0.1
                            p0 = [A_0,mu_0,sigma_0] #C_0]
                            #changed from .03 to .04 on 10-02-19
                            c1 = .04#c1 control the constraint on the line center
                            #limit the peak height unless the peak is saturated
                            if (data_in[frame,location] == 65535):
                                count_wiggle = 200000
                            else:
                                count_wiggle = 100                           
                            bounds_init = ((data_in[frame,location] - 100 , peak_ls[1][frame,j,0] - c1,0), \
                                           (data_in[frame,location] + count_wiggle, peak_ls[1][frame,j,0] + c1, 0.14))
                            #fit to a gaussian using the curve_fit module:
                            coeff, var_matrix = curve_fit(gaussian, data_in[0,:], data_in[frame,:], p0=p0, \
                                bounds = bounds_init) 
                            fit = gaussian(data_in[0,:],*coeff)
                            #plt.plot(data_in[0,:], fit)
                            #%%
                            #intersections between gaussian and baseline                    
                            intersections = np.argwhere(np.diff(np.sign(fit - baseline))).flatten()
                            #new x limits are the intersections between baseline and gaussian functions
                            lower_x = data_in[0,intersections[0]]
                            upper_x = data_in[0,intersections[1]]
                            #background counts for this peak:
                            bkg_counts = quad(poly,lower_x,upper_x)
                            #numerically integrated intens:
                            int_intens = quad(gaussian, lower_x , upper_x, args = (coeff[0],coeff[1],coeff[2]))
                            #analytically integrated intensity; save both later
                            analytic_intens = analytic_counts(lower_x, upper_x, coeff[0],coeff[1],coeff[2])
                            #%linspace integration; outdated but keeping for future if needed]
                            #x_data = data_in[0,:]
                            #x_range = np.linspace(data_in[0,location] - 0.5, data_in[0,location] + 0.5 ,num=10000, dtype = float)
                            #linspace_intens = quad(gaussian, x_range[0] , max(x_range), args = (coeff[0],coeff[1],coeff[2]))
                            #print('Numerical Intens: {:}'.format(int_intens[0] - bkg_counts[0]))
                            #print('Analytic Intens: {:}'.format(analytic_intens - bkg_counts[0]))
                            #print(coeff[1])
                            #print(coeff[2])
                            #%find the total gaussian fit function & calc. the synth. spectra
                            total_gauss = gauss_sum(total_gauss, data_in[0,:], coeff[0],coeff[1],coeff[2])
                            synthetic_spec = synthetic_spectra(synthetic_spec, fit[:],intersections)
                            #save wavelength, numerical intens, analytic intens, and frame # for each peak
                            temp[0,0] = lines_in_region[i]
                            #%%
                            temp[0,1] = int_intens[0] - bkg_counts[0]
                            if (temp[0,1] < 0):
                                temp[0,1] = 0
                            temp[0,2] = analytic_intens - bkg_counts[0]
                            temp[0,3] = frame 
                            intens_arr = np.append(intens_arr,temp,axis=0)
                        except RuntimeError:
                            print('Error - curve_fit failed for peak at {:}'.format(lines_in_region[i]))
                           #file.close()
                        except IndexError:
                            print('Error - No Intersections (peak {:} DNE)'.format(round(lines_in_region[i],5)))   
        #Commented out on 10-02-19 as a precautionary measure. Will re-enable if necessary
        """
            if (lines_in_region not in temp):
                #try to fit a peak in case the peak was not detected. If it works, save it, otherwise trash it
                try:
                    #initial guesses for the fit parameters:
                    location = find_nearest_index(data_in[0,:], lines_in_region[i])
                    A_0 = data_in[frame,location]
                    mu_0 = data_in[0,location]
                    sigma_0 = 0.1
                    p0 = [A_0,mu_0,sigma_0] #C_0]
                    c1 = .02#c1 control the constraint on the line center
                    count_wiggle = 100 
                    #wiggle room in the counts = count_wiggle. Creative variable name
                    bounds_init = ((data_in[frame,location] - 100 , data_in[0,location] - c1,0), \
                                           (data_in[frame,location] + count_wiggle, data_in[0,location] + c1, 0.14))                          
                    coeff, var_matrix = curve_fit(gaussian, data_in[0,:], data_in[frame,:], p0=p0, \
                                bounds = bounds_init) 
                    fit = gaussian(data_in[0,:],*coeff)       
                    #find where the baseline function & gaussian intersect:
                    intersections = np.argwhere(np.diff(np.sign(fit - baseline))).flatten()
                    #save the upper and lower bounds of the region containing the peak above the basline
                    lower_x = data_in[0,intersections[0]]
                    upper_x = data_in[0,intersections[1]]
                    #numerically integrate the background profile:
                    bkg_counts = quad(poly,lower_x,upper_x)
                    int_intens = quad(gaussian, lower_x , upper_x, args = (coeff[0],coeff[1],coeff[2]))
                    analytic_intens = analytic_counts(lower_x, upper_x, coeff[0],coeff[1],coeff[2])
                    total_gauss = gauss_sum(total_gauss, data_in[0,:], coeff[0],coeff[1],coeff[2])
                    synthetic_spec = synthetic_spectra(synthetic_spec, fit[:],intersections)
                    #save wavelength, numerical intens, analytic intens, and frame # for each peak
                    temp[0,0] = lines_in_region[i]
                    #save the integrated intensity minus the background contribution:
                    temp[0,1] = int_intens[0] - bkg_counts[0]
                    temp[0,2] = analytic_intens - bkg_counts[0]
                    temp[0,3] = frame 
                    intens_arr = np.append(intens_arr,temp,axis=0)
                except RuntimeError:
                    print('Error - curve_fit failed for peak at {:}'.format(lines_in_region[i]))
                except IndexError:
                    print('Error - No Intersections (peak {:} DNE)'.format(round(lines_in_region[i],5)))
        """
    intensity_out = np.empty((0,2))
    os.chdir(cwd)
    #sort the output and return the integrated intensities
    for i in range(0,len(lines[:])):
        temp_intens = np.empty((1,2))
        temp_intens[0,0] = lines[i]
        for j in range(0,len(intens_arr[:,0])):
            if (intens_arr[j,0] == lines[i]):
                temp_intens[0,1] = temp_intens[0,1] + intens_arr[j,1]
        intensity_out = np.append(intensity_out, temp_intens, axis=0)
    #SJB 09-06-19 I know the 'True' and 'absolute' functions work. 
    if ((normalize == 'True') and (method == 'absolute')):
        max_val = np.amax(intensity_out[:,1])
        intensity_out[:,1] = scale_max * intensity_out[:,1] / max_val  
    #work with a copy of intensity_out. We'll use this copy to change the values in intensity_out to their normalized values
    #later
    intensity_out_temp = intensity_out
    #Not working
    if ((normalize == 'True') and (method == 'relative')):
        for r in range(0,len(regions[:,0])):
            lower_x = regions[r,1]
            upper_x = regions[r,2]
            region_list = np.empty((0,2))
            #find the lines within each wavelength region:
            for i in range(0,len(intensity_out_temp[:,0])):
                temp = np.empty((1,2))
                if (lower_x < intensity_out_temp[i,0] < upper_x):
                    temp[0,0] = intensity_out_temp[i,0]
                    temp[0,1] = intensity_out_temp[i,1]
                    region_list = np.append(region_list,temp, axis = 0)
            
            if (r == 10):
                np.savetxt('region_list_test.txt', region_list)
            #find the maximum for lines in this wavelength region 
            try:
                max_val = np.amax(region_list[:,1])
                #normalize to the maximum value in this region
                region_list[:,1] = scale_max * region_list[:,1] / max_val
               
                #replace the intensities in the master list with the normalized values:
                for i in range(0,len(intensity_out_temp[:,0])):
                        j_max = len(region_list)    
                        for j in range(0,len(region_list)):
                            #Added on 10-15-19 to fix issue with blend intensity                                
                            if (region_list[j,0] == intensity_out_temp[i,0]):
                                if (len(region_list) < 2):
                                    intensity_out_temp[i,1] = region_list[j,1]
                                if (region_list[j,0] != region_list[j-1,0]):
                                    intensity_out_temp[i,1] = region_list[j,1]

            except ValueError:
                #STEVE I should add in check that if the peak IS NOT found then it adds a 0 to the array
                #will do this next time so the entire intensity code is self-contained
                #and won't rely on any post processing or checking
                print('No Peaks in wavelength range {:} : {:} - {:} nm'.format(countr,regions[r,1],regions[r,2]))
    intensity_out = intensity_out_temp
    return intensity_out



def intensity_search(array,frame1,frame2,normalize,scale_max):
    ########################################################################################
    ######################  Explanation of the parameters   ################################
    #                                                                                      # 
    #   array: input array of line positions                                               # 
    #   frame1: start of frame range over which to search                                  # 
    #   frame2: end of frame range over which to search                                    #  
    #   normalize: boolean; if 'True', normalize the largest value to scale_max            #
    #   scale_max: absolute maximum by which the found line intensities are normalized     #
    #                                                                                      #
    ########################################################################################
    cwd = os.getcwd()
    #define the data directory (probe2), output dir. for plots (frame plots), and shot_list (wav. regions)
    probe2 = 'C:/Users/steve_000/Desktop/Research_Files/PRA_gold_spectroscopy/raw_data_noxrays_nosum/2nd_probe'
    frame_plots = 'C:/Users/steve_000/Desktop/Research_Files/PRA_gold_spectroscopy/07-15-19/frame_test_fits'
    #LOPT_folder = 'C:/Users/steve_000/Desktop/Research_Files/PRA_gold_spectroscopy/07-15-19/Au_I_LOPT'
    shot_list = genfromtxt('shots.csv',delimiter = ',')
    #lines file to look for; pull from LOPT folder:
    #os.chdir(LOPT_folder)
    #au_I_lines_dat = glob.glob('au_I_lines*.txt')
    #lines_in = np.genfromtxt(au_I_lines_dat[0],delimiter ='\t', dtype = str,skip_header=1)#final output array:
    lines = array
    os.chdir(cwd)
    temp = np.empty((1))
    #initial output of intensity search is  intensity is intens_arr
    intens_arr = np.empty((0,4))
    for countr in range(0,17):
        #load the data for the shot in wavelength range defined by countr:
        os.chdir(probe2)
        probe2_6cm = str(int(shot_list[countr,4])) + '.csv'
        probe2_6cm_dat = genfromtxt(probe2_6cm,delimiter =',')
        probe2_6cm_shot = str(int(shot_list[countr,4]))
        probe2_6cm_num = int(probe2_6cm_shot)
        data_in = probe2_6cm_dat
        x_min = data_in[0,0]
        x_max = data_in[0,1339]
        lines_in_region = np.empty((0))
        #temp arr for holding lines
        temp1 = np.empty((1))
        #grab the lines within this wavelength window
#
        #%%
        for i in range(0,len(lines)):    
            if (x_min < lines[i] < x_max):
                temp1 = lines[i]
                lines_in_region = np.append(lines_in_region,temp1)
        lines_in_region = lines_in_region[lines_in_region[:].argsort()]                
        #%%
        #%Find the peaks using my standard peak finding function:
        #frame1 = 60
        #frame2 = 120
        #peak list from peak finding code. Do for all frames from frame1 -frame2
        peak_ls = peak_find(frame1,frame2,data_in)
        os.chdir(frame_plots)
        #loop over all frames :)
        for frame in range(frame1, frame2):
            #Find the baseline values & initialize values for the total gaussian and the spectra function
            total_gauss = 0
            #modified the peakutils baseline function to return the coefficients too :)
            deg = 10
            baseline_tot = baseline_mod(data_in[frame,:], deg, max_it=100, tol=1e-5)
            baseline = baseline_tot[0]
            #define an array of coefficients for poly1d function:
            #fit our x values to a 'deg' degree polynomial for the baseline y values; return coeffs
            poly_params = np.polyfit(data_in[0,:], baseline, deg)
            #define a polynomial with the resulting coeffs. Use the function 'poly' below for integrating
            poly = np.poly1d(poly_params)
            poly_dat = poly(data_in[0,:])
            #the baseline acts as the starting point for constructing the 'synthetic spectra'
            synthetic_spec = baseline
            for i in range(0,len(lines_in_region)):
                #open up an output file for troubleshooting if need be:
                #file = open('testfile_frame_{:}.txt'.format(lines_in_region), 'w')
                #tempa array to hold each peaks data before saving:
                temp = np.empty((1,4))
                #only look for the peak & fit IF the peak is in this wavelength region:
                location = find_nearest_index(data_in[0,:], lines_in_region[i])
                for j in range(0,len(data_in[0,:])):
                    #check to make sure the peak finding found the peak in this frame:
                    if (abs(lines_in_region[i] - peak_ls[1][frame,j,0]) < .04):
                        #make sure the peak isn't at the borders:
                        if (5 < location < 1334):
                            #try to fit; otherwise throw an error and go to next peak:
                            try:
                                #Set the initial parameters for height, center, and width
                               #%%
                                #location = 1057
                                A_0 = data_in[frame,location]
                                #mu_0 = data_in[0,location]
                                mu_0 = peak_ls[1][frame,j,0]
                                sigma_0 = 0.1
                                p0 = [A_0,mu_0,sigma_0] #C_0]
                                c1 = .03#c1 control the constraint on the line center
                                #limit the peak height unless the peak is saturated
                                if (data_in[frame,location] == 65535):
                                    count_wiggle = 200000
                                else:
                                    count_wiggle = 100                           
                                bounds_init = ((data_in[frame,location] - 100 , peak_ls[1][frame,j,0] - c1,0), \
                                               (data_in[frame,location] + count_wiggle, peak_ls[1][frame,j,0] + c1, 0.14))
                                #fit to a gaussian using the curve_fit module:
                                coeff, var_matrix = curve_fit(gaussian, data_in[0,:], data_in[frame,:], p0=p0, \
                                    bounds = bounds_init) 
                                fit = gaussian(data_in[0,:],*coeff)
                                #plt.plot(data_in[0,:], fit)
                                #%%
                                #intersections between gaussian and baseline                    
                                intersections = np.argwhere(np.diff(np.sign(fit - baseline))).flatten()
                                #new x limits are the intersections between baseline and gaussian functions
                                lower_x = data_in[0,intersections[0]]
                                upper_x = data_in[0,intersections[1]]
                                #background counts for this peak:
                                bkg_counts = quad(poly,lower_x,upper_x)
                                #numerically integrated intens:
                                int_intens = quad(gaussian, lower_x , upper_x, args = (coeff[0],coeff[1],coeff[2]))
                                #analytically integrated intensity; save both later
                                analytic_intens = analytic_counts(lower_x, upper_x, coeff[0],coeff[1],coeff[2])
                                #%linspace integration; outdated but keeping for future if needed]
                                #x_data = data_in[0,:]
                                #x_range = np.linspace(data_in[0,location] - 0.5, data_in[0,location] + 0.5 ,num=10000, dtype = float)
                                #linspace_intens = quad(gaussian, x_range[0] , max(x_range), args = (coeff[0],coeff[1],coeff[2]))
                                #print('Numerical Intens: {:}'.format(int_intens[0] - bkg_counts[0]))
                                #print('Analytic Intens: {:}'.format(analytic_intens - bkg_counts[0]))
                                #print(coeff[1])
                                #print(coeff[2])
                                #%find the total gaussian fit function & calc. the synth. spectra
                                total_gauss = gauss_sum(total_gauss, data_in[0,:], coeff[0],coeff[1],coeff[2])
                                synthetic_spec = synthetic_spectra(synthetic_spec, fit[:],intersections)
                                #save wavelength, numerical intens, analytic intens, and frame # for each peak
                                temp[0,0] = lines_in_region[i]
                                temp[0,1] = int_intens[0] - bkg_counts[0]
                                temp[0,2] = analytic_intens - bkg_counts[0]
                                temp[0,3] = frame 
                                intens_arr = np.append(intens_arr,temp,axis=0)
                            except RuntimeError:
                                print('Error - curve_fit failed for peak at {:}'.format(lines_in_region[i]))
                               #file.close()
                            except IndexError:
                                print('Error - No Intersections (peak {:} DNE)'.format(round(lines_in_region[i],5)))   
                if (lines_in_region not in temp):
                    #try to fit a peak in case the peak was not detected. If it works, save it, otherwise trash it
                    try:
                        location = find_nearest_index(data_in[0,:], lines_in_region[i])
                        A_0 = data_in[frame,location]
                        mu_0 = data_in[0,location]
                        sigma_0 = 0.1
                        p0 = [A_0,mu_0,sigma_0] #C_0]
                        c1 = .02#c1 control the constraint on the line center
                        count_wiggle = 100 
                        bounds_init = ((data_in[frame,location] - 100 , data_in[0,location] - c1,0), \
                                               (data_in[frame,location] + count_wiggle, data_in[0,location] + c1, 0.14))                          
                        coeff, var_matrix = curve_fit(gaussian, data_in[0,:], data_in[frame,:], p0=p0, \
                                    bounds = bounds_init) 
                        fit = gaussian(data_in[0,:],*coeff)       
                        intersections = np.argwhere(np.diff(np.sign(fit - baseline))).flatten()
                        lower_x = data_in[0,intersections[0]]
                        upper_x = data_in[0,intersections[1]]
                        bkg_counts = quad(poly,lower_x,upper_x)
                        int_intens = quad(gaussian, lower_x , upper_x, args = (coeff[0],coeff[1],coeff[2]))
                        analytic_intens = analytic_counts(lower_x, upper_x, coeff[0],coeff[1],coeff[2])
                        total_gauss = gauss_sum(total_gauss, data_in[0,:], coeff[0],coeff[1],coeff[2])
                        synthetic_spec = synthetic_spectra(synthetic_spec, fit[:],intersections)
                        #save wavelength, numerical intens, analytic intens, and frame # for each peak
                        temp[0,0] = lines_in_region[i]
                        temp[0,1] = int_intens[0] - bkg_counts[0]
                        temp[0,2] = analytic_intens - bkg_counts[0]
                        temp[0,3] = frame 
                        intens_arr = np.append(intens_arr,temp,axis=0)
                    except RuntimeError:
                        print('Error - curve_fit failed for peak at {:}'.format(lines_in_region[i]))
                    except IndexError:
                        print('Error - No Intersections (peak {:} DNE)'.format(round(lines_in_region[i],5)))
    intensity_out = np.empty((0,2))
    os.chdir(cwd)
    #sort the output and return the integrated intensities
    for i in range(0,len(lines[:])):
        temp_intens = np.empty((1,2))
        temp_intens[0,0] = lines[i]
        for j in range(0,len(intens_arr[:,0])):
            if (intens_arr[j,0] == lines[i]):
                temp_intens[0,1] = temp_intens[0,1] + intens_arr[j,1]
        intensity_out = np.append(intensity_out, temp_intens, axis=0)
    
    if (normalize == 'True'):
        max_val = np.amax(intensity_out[:,1])
        intensity_out[:,1] = scale_max * intensity_out[:,1] / max_val            
            
    return intensity_out

#%%

##########

#Added pon 07-30-19
def baseline_mod(y, deg=None, max_it=None, tol=None):
    """
    Computes the baseline of a given data.

    Iteratively performs a polynomial fitting in the data to detect its
    baseline. At every iteration, the fitting weights on the regions with
    peaks are reduced to identify the baseline only.
    This is a modified version of the PeakUtils baseline
    """
    # for not repeating ourselves in `envelope`
    if deg is None: deg = 3
    if max_it is None: max_it = 100
    if tol is None: tol = 1e-3
    
    order = deg + 1
    coeffs = np.ones(order)

    # try to avoid numerical issues
    cond = math.pow(abs(y).max(), 1. / order)
    x = np.linspace(0., cond, y.size)
    base = y.copy()

    vander = np.vander(x, order)
    vander_pinv = LA.pinv2(vander)

    for _ in range(max_it):
        coeffs_new = np.dot(vander_pinv, y)

        if LA.norm(coeffs_new - coeffs) / LA.norm(coeffs) < tol:
            break

        coeffs = coeffs_new
        base = np.dot(vander, coeffs)
        y = np.minimum(y, base)

    return base, coeffs



##Added on 07-30-19 to steve.py
def synthetic_spectra(spec, fit, intersections):
    #note: to use this, define the synthetic_seed as the baseline function at the beginning of the run
    #generates a synthetic spectra of baseline + gaussians over top
    for i in range(intersections[0] + 1, intersections[1] + 1):
        spec[i] = fit[i]
    return spec   

#Note: to use gauss_sum function, a "seed = 0" must be used to start the sum when looping over peaks
def gauss_sum(seed,x,a,mu,sigma):
    seed = seed + a*np.exp(-(x-mu)**2 / (2*sigma**2))
    return seed

#standard gaussian function:
def gaussian(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2*sigma**2))

#Returns analytic solution to integral of a gaussian between two limits given by lower_x and upper_x
def analytic_counts(lower_x, upper_x, A, mu, sigma):
    lower = -(A * sigma) * sqrt(np.pi/2) * erf((1/sigma) * (mu - lower_x) / sqrt(2))
    upper = -(A * sigma) * sqrt(np.pi/2) * erf((1/sigma) * (mu - upper_x) / sqrt(2))
    diff = upper - lower
    return diff


#############################################
#Added 07-03-19, SJB
#def wavenumber_to_nm(val):
#    nm = 1e9 * 1.98644568e-25 / (val * 1.9863E-23)
 #   return(nm)

def rel_intens(array,frame1,frame2,shot_list):
    cwd = os.getcwd()
    probe2 = 'C:/Users/steve_000/Desktop/Research_Files/PRA_gold_spectroscopy/raw_data_noxrays_nosum/2nd_probe'
    peak_list_initial = array
    #shot_list = genfromtxt('shots.csv',delimiter = ',')
    scatter_data_out = np.empty((0,4))
    for countr in range(0,17):
        os.chdir(probe2)
        probe2_6cm = str(int(shot_list[countr,4])) + '.csv'
        probe2_6cm_dat = genfromtxt(probe2_6cm,delimiter =',')
        min_lambda = probe2_6cm_dat[0,0]
        max_lambda = probe2_6cm_dat[0,1339]
        peak_list = np.empty((0))
          
        for i in range(0,len(peak_list_initial)):
            if (min_lambda < peak_list_initial[i] < max_lambda):
                peak_list = np.append(peak_list,peak_list_initial[i])         
        peak_list = peak_list[np.argsort(peak_list)]
        #peak find:
        pk_find_probe2_6cm = peak_find(frame1,frame2,probe2_6cm_dat)    
        #scatter data:
        for i in range(0,len(peak_list)):
            scatter_probe2_6cm = scatter_data(probe2_6cm_dat,pk_find_probe2_6cm[1],peak_list[i],frame1,frame2)
            scatter_data_out = np.append(scatter_data_out,scatter_probe2_6cm,axis=0)
           
    os.chdir(cwd)
    for i in range(0,len(scatter_data_out[:,0])):
        scatter_data_out[i,3] = scatter_data_out[i,2] / 65.535
    
    addon = np.empty((0,4))
    for i in range(0,len(peak_list_initial)):
        presence = 1
        temp = np.zeros((1,4))
        for j in range(0,len(scatter_data_out[:,0])):    
            if (peak_list_initial[i] == scatter_data_out[j,0]):
                presence = 0
        if (presence == 1):
            temp[0,0] = peak_list_initial[i]
            temp[0,1] = 0
            temp[0,2] = 0
            temp[0,3] = 1
            addon = np.append(addon,temp,axis=0)
    scatter_data_out = np.append(scatter_data_out,addon,axis=0)
    scatter_data_out = scatter_data_out[scatter_data_out[:,0].argsort()]    
        
    array2 = scatter_data_out     
    return(array2)



#Added: ??
def read_prin(fil):
    f = open(fil)
    data = np.loadtxt(f,delimiter=',')
    inten = np.reshape(data[:,1],(-1,1340))
    return np.vstack([data[0:1340,0],inten])

#To save space in peak_finding codes
def stdev(y):
   return(np.std(y))

#Adding stdevs in quadrature; used in peak finding binning procedure
#Added: ??
def quad_std(x):
    sum = 0
    for i in range(0,len(x)):
        sum = sum + (x[i])**2
    mean = np.sqrt(sum)
    return(mean)  

#Added: 03-27-19    
def find_nearest(array, value):
    idx = (np.abs(array - value)).argmin()
    return idx

def find_nearest_index(array,value):
    diff_temp = abs(value-array[0])
    index = 0 
    for i in range(0,len(array[:])):
        diff = abs(array[i] - value)
        if (diff < diff_temp):
            index = i
            diff_temp = diff
    return index

#Used to generate the data required to make a position scatter plot
def scatter_data(raw_data,all_peaks_out,peak,frame1,frame2):
    #Note to self, restructure this to 3d to store all peak information?
    peak_intens = np.empty((0,4))       
    #UP TO HERE WORKS
    for frame in range(frame1,frame2):
        #find nearest peak to the peak given in peak list and compare to each interp_peak
        frame_peak = find_nearest(all_peaks_out[frame,:,0],peak)
        #find index of that peak in the data_in array
        peak_index = find_nearest(raw_data[0,:], all_peaks_out[frame,frame_peak,0])
        if (abs(raw_data[0,peak_index] - peak) < .05):
            #If the nearest peak is within some threshold, i.e. 0.04nm, then
            peak_temp = np.empty((1,4))
            peak_temp[0,0] = peak #initial peak
            peak_temp[0,1] = all_peaks_out[frame,frame_peak,0] #interp. value of peak for record keeping
            peak_temp[0,2] = all_peaks_out[frame,frame_peak,1] #counts
            peak_temp[0,3] = frame
            peak_intens = np.append(peak_intens,peak_temp,axis=0) 
    return(peak_intens)

#Peak finding function; in the future, I should modify this to be user-adjustable from the function call
#Added: 03-27-19
def peak_find(frame1,frame2,array):
    array2 = np.zeros((121,1340,2))
    for frame in range(frame1,frame2):#range(frame1,frame2):
        min_lambda = array[0,0]
        max_lambda = array[0,1339]
        baseline = peakutils.baseline(array[frame,:], deg=3, max_it=100, tol=1e-5)
        indexes = peakutils.peak.indexes(array[frame,:], thres=baseline + 100, min_dist=2, thres_abs=True)
        interp_peaks = peakutils.peak.interpolate(array[0,:],array[frame,:],indexes,width = 1)
        numPeaks = len(interp_peaks[:])
        array2[0,frame,0] = frame
        k=0
        for j in range(0,numPeaks):
            if (min_lambda < interp_peaks[j] < max_lambda):
                array2[frame,k,0] = interp_peaks[j]
                array2[frame,k,1] = array[frame,indexes[j]]
                k = k +1
    return(array,array2)

#Simple vac to air conversion. Taken from a ~2002 ApJ paper (I think)
#Note: this is also the function colradpy uses :)
#Added: 03-27-19
def vac_to_air(array):
    #Added on 03-05-19
    #Edited on 04-11-19 to only require input of a single list
    #NOTE: inputs must be in angstroms
    for i in range(0,len(array)):
        s2 = 1e4 / array[i]
        n = 1 + 0.0000834254 + (0.02406147)/(130 - s2) + (0.00015998)/(38.9 - s2)
        array[i] = array[i] / n

    return(array)
    
def val_vac_to_air(val):
    #Added on 03-05-19
    #Edited on 04-11-19 to only require input of a single list
    #NOTE: inputs must be in angstroms
    s2 = 1e4 / val
    n = 1 + 0.0000834254 + (0.02406147)/(130 - s2) + (0.00015998)/(38.9 - s2)
    val = val / n

    return(val)    
    
#??
def trim_by_column(column,value,array):
    logic = 0
    i = 0 
    while (logic < 10):
        if (i == len(array[:,0])):
            break
        if (array[i,column] < value):
            array = np.delete(array,(i),axis=0)
            i = i - 1
        i = i + 1
    return(array)

def trim_by_counts(counts,array):
    logic = 0
    i = 0 
    while (logic < 10):
        if (i == len(array[:,0])):
            break
        if (array[i,1] < counts):
            array = np.delete(array,(i),axis=0)
            i = i - 1
        i = i + 1
    return(array)
    
def trim_by_number_peaks(peaks,array):
    logic = 0
    i = 0 
    while (logic < 10):
        if (i == len(array[:,0])):
            break
        if (array[i,2] < peaks):
            array = np.delete(array,(i),axis=0)
            i = i - 1
        i = i + 1
    return(array)
    
#If save_all = 1, saves final line list at end; if save_each =, saves individual line lists with all IDs
def NIST_comparison(array,NIST_path,output_path,thresh,save_all,save_each):
    temp_dir = os.getcwd()
    os.chdir(NIST_path)
    NIST_files = glob.glob('*.txt')
    os.chdir(temp_dir)
    IDs = np.empty(len(array[:,0]),dtype ='<U256')
    start=time.time()
    for i in range(0,len(IDs)):
        IDs[i] = ''
    temp_species = np.empty((1,len(array[0,:])))
    for file in NIST_files:
        os.chdir(NIST_path)
        data_species = np.empty((0,len(array[0,:])))
        lines_file = open(file)  
        lines = np.loadtxt(lines_file)
        lines_ID = file.rsplit('.',1)[0]
        for i in range(0,len(array[:,0])):   
            #for each PEAK detected, check against ALL NIST lines:
            for k in range(0,len(lines[:])):
                #Add 'tag' to ID array IF the NIST line is within +/- thresh from our peaks
                if (abs(array[i,0] - lines[k]) < thresh): 
                    IDs[i] = IDs[i] + lines_ID + ','
                    for n in range (0,len(array[0,:])):
                        temp_species[0,n] = array[i,n]
                    data_species = np.append(data_species,temp_species,axis=0)
        if (save_each == 'true'):
            os.chdir(output_path)
            np.savetxt('{:}.txt'.format(lines_ID),data_species,delimiter = ',')
            os.chdir(NIST_path)
        lines_file.close()
#going to need to generalize this at some point
    a = array[:,0]
    b = array[:,1]
    c = array[:,2]
    d = array[:,3]
    e = array[:,4]
    f2 = array[:,5]
    g = array[:,6]
    h = array[:,7]
    i = array[:,8]
    j = array[:,9]
    k = array[:,10]
    l = IDs[:]
    if (save_all == 'true'):
        os.chdir(output_path)
        with open('line_comp_out with ID.csv','w') as f:
                writer = csv.writer(f, delimiter=',')
                writer.writerow(["Lambda", "Stdev", "# Counts", "# Peaks", "Au Out","Ni Out","Au 3cm", "Ni 3cm", "Au 6cm", "Ni 6cm", "Au 9cm", "NIST ID"])
                writer.writerow(["See 'Classifications' file for comparison meanings"])
                writer.writerows(zip(a,b,c,d,e,f2,g,h,i,j,k,l))   
        temp_new_lines = np.empty((1,len(array[0,:])))
        new_lines = np.empty((0,len(array[0,:])))
        for i in range(0,len(array[:,0])):
            if (IDs[i] == ''):
                for j in range(0,len(temp_new_lines[0,:])):
                    temp_new_lines[0,j] = array[i,j]
                new_lines = np.append(new_lines,temp_new_lines,axis=0)        
        a = new_lines[:,0]
        b = new_lines[:,1]
        c = new_lines[:,2]
        d = new_lines[:,3]
        e = new_lines[:,4]
        f2 = new_lines[:,5]
        g = new_lines[:,6]
        h = new_lines[:,7]
        i = new_lines[:,8]
        j = new_lines[:,9]
        k = new_lines[:,10]
        with open('New Lines with no ID.csv','w') as f:
                writer = csv.writer(f, delimiter=',')
                writer.writerow(["Lambda", "Stdev", "# Counts", "# Peaks", "Au Out","Ni Out","Au 3cm", "Ni 3cm", "Au 6cm", "Ni 6cm", "Au 9cm", "NIST ID"])
                writer.writerow(["See 'Classifications' file for comparison meanings"])
                writer.writerows(zip(a,b,c,d,e,f2,g,h,i,j,k))   
    
                
        
    os.chdir(temp_dir)
    print ('First list, all lines with ID: {:} s'.format(time.time()-start))
    #Analyze each element + charge individually:
    if (save_each == 'true'):
        os.chdir(output_path)
        species_files = glob.glob('*.txt')
        for input in species_files:
          start=time.time()
          data_in = genfromtxt(input,delimiter = ',')
          IDs2 = np.empty(len(array[:,0]),dtype ='<U256')
          input_ID = input.rsplit('.',1)[0]
          for i in range(0,len(IDs2)):
              IDs2[i] = ''
          os.chdir(NIST_path)
          for file in NIST_files:
              lines_file = open(file)  
              lines = np.loadtxt(lines_file)
              lines_ID = file.rsplit('.',1)[0]
              for i in range(0,len(data_in[:,0])):   
                  for k in range(0,len(lines[:])):
                      if (abs(data_in[i,0] - lines[k]) < thresh): 
                          IDs2[i] = IDs2[i] + lines_ID + ','
          os.chdir(output_path)
          a = data_in[:,0]
          b = data_in[:,1]
          c = data_in[:,2]
          d = data_in[:,3]
          e = data_in[:,4]
          f2 = data_in[:,5]
          g = data_in[:,6]
          h = data_in[:,7]
          i = data_in[:,8]
          j = data_in[:,9]
          k = data_in[:,10]
          l = IDs2[:]
          with open('{:} with all IDs.csv'.format(input_ID),'w') as f:
              writer = csv.writer(f, delimiter=',')
              writer.writerow(["Lambda", "Stdev", "# Counts", "# Peaks", "Au Out","Ni Out","Au 3cm", "Ni 3cm", "Au 6cm", "Ni 6cm", "Au 9cm", "NIST ID"])
              writer.writerow(["See 'Classifications' file for comparison meanings"])
              writer.writerows(zip(a,b,c,d,e,f2,g,h,i,j,k,l))    
                  
          print ('{:} : {:} seconds'.format(input_ID,time.time()-start))
    os.chdir(temp_dir)
    return()  
    
def line_ratio_comparison(array1,array2,thresh):
    lr_array = np.zeros((len(array1[:,0]),1))
    for i in range(0,len(array1[:,0])):
        presence = 0
        avg_ar1 = (array1[i,1] / array1[i,2])
        for j in range(0,len(array2[:,0])):
            avg_ar2 = (array2[j,1] / array2[j,2])
            lr = avg_ar1 / avg_ar2
            #Compare avg counts (counts / # peaks) between 2 arrays
            if (abs(array1[i,0] - array2[j,0]) < thresh):
                if (abs(lr) < 1):
                    lr_array[i] = 1
                    presence = 1 #weaker than probe out data
                if (1 < abs(lr) < 2):  
                    lr_array[i] = 2 #comparable to probe out data
                    presence = 1
                if (2 < abs(lr) < 5): #greater than probe out data
                    lr_array[i] = 3 
                    presence = 1
                if (abs(lr) > 5):
                    lr_array[i] = 4 #much stronger than probe out data
                    presence = 1    
        if (presence == 0):
            lr_array[i] = 0
    #SB: Should I add this option?
    #Option to save this step of the process
    #np.savetxt('Total Line List with IDs.txt',data_total,fmt='%.7f',delimiter=',')
    return(lr_array[:])