"""
SJB
Peak Intensity Code
"""
import os
import numpy as np
import scipy
from scipy import stats
from scipy.optimize import curve_fit
from numpy import inf
import peakutils
import csv
import glob
import re
from numpy import genfromtxt
from scipy.integrate import *
from math import erf
from math import *
################################################
#Peak Finding Function
def peak_find(frame1,frame2,array):
    array2 = np.zeros((121,1340,2))
    for frame in range(frame1,frame2):#range(frame1,frame2):
        min_lambda = array[0,0]
        max_lambda = array[0,1339]
        #Estimate background counts
        baseline = peakutils.baseline(array[frame,:], deg=3, max_it=100, tol=1e-5)
        #Approximate peak positions
        indexes = peakutils.peak.indexes(array[frame,:], thres=baseline + 100,\
                                         min_dist=2, thres_abs=True)
        #Refine peak positions
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

def intensity_search_v3(array,frame_in,normalize,scale_max,method,regions, \
                        data_directory,frame_source,calib,correct):
    ##########################################################################################
    ######################  Explanation of the parameters   ##################################
    #                                                                                                                                          
    #   array: input array of line positions                                                                                                   
    #   frame: frame you want the intensities from                                                                                             
    #   normalize: boolean; if 'True', normalize the largest value to scale_max                                                               
    #   scale_max: absolute maximum by which the found line intensities are normalized                                                        
    #   method: select either 'relative' or 'absolute'                                                                                         
    #         'relative': normalizes to scale_max within each wavelength region                                                               
    #         'absolute': normalizes to scale_max over all wavelength ranges (not reliable)                                                   
    #   regions: array containing the wavelength regions from the raw data. Ex:                                                               
    #            Format: (region #) (lower lambda) (upper lambda) (frame for intensities) 
    #            (deg. for baseline polynomial)                      
    #   frame_source: 'user' uses the value of 'frame' in the function call, or 'file' 
    #    uses the col #4 from regions array                     
    #   calib: intensity response of detector on 0 - 100 scale; 2 col format: 
    #   Wavelength (nm) | Response. (defined 0 as first element)        #
    #   correct: boolean controlling whether or not to re-calibration with the above calib file                                               
    ###########################################################################################
    cwd = os.getcwd()
    #define the data directory (probe2), output dir. for plots (frame plots), and 
    #shot_list (wav. regions)
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
        #Find the baseline values & initialize values for the total gaussian and 
        #the spectra function
        total_gauss = 0
        #modified the peakutils baseline function to return the coefficients too :)
        baseline_tot = baseline_mod(data_in[frame,:], deg, max_it=100, tol=1e-5)
        baseline = baseline_tot[0]
        #define an array of coefficients for poly1d function:
        #fit our x values to a 'deg' degree polynomial for the baseline y values; return coeffs
        poly_params = np.polyfit(data_in[0,:], baseline, deg)
        #define a polynomial with the resulting coeffs. Use the function 'poly' below 
        #for integrating
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
                            bounds_init = ((data_in[frame,location] - 100 , \
                                            peak_ls[1][frame,j,0] - c1,0), \
                                           (data_in[frame,location] + count_wiggle, \ 
                                            peak_ls[1][frame,j,0] + c1, 0.14))
                            #fit to a gaussian using the curve_fit module:
                            coeff, var_matrix = curve_fit(gaussian, data_in[0,:], \ 
                                data_in[frame,:], p0=p0, \
                                bounds = bounds_init) 
                            fit = gaussian(data_in[0,:],*coeff)
                            #intersections between gaussian and baseline                    
                            intersections = \
							np.argwhere(np.diff(np.sign(fit - baseline))).flatten()
                            #new x limits are the intersections between baseline 
							#and gaussian functions
                            lower_x = data_in[0,intersections[0]]
                            upper_x = data_in[0,intersections[1]]
                            #background counts for this peak:
                            bkg_counts = quad(poly,lower_x,upper_x)
                            #numerically integrated intens:
                            int_intens = quad(gaussian, lower_x , upper_x, \
							args = (coeff[0],coeff[1],coeff[2]))
                            #analytically integrated intensity; save both later
                            analytic_intens = \
							analytic_counts(lower_x, upper_x, coeff[0],coeff[1],coeff[2])
                            #save wavelength, numerical intens, analytic intens, 
							#and frame # for each peak
                            temp[0,0] = lines_in_region[i]
                            temp[0,1] = int_intens[0] - bkg_counts[0]
                            if (temp[0,1] < 0):
                                temp[0,1] = 0
                            temp[0,2] = analytic_intens - bkg_counts[0]
                            temp[0,3] = frame 
                            intens_arr = np.append(intens_arr,temp,axis=0)
                        except RuntimeError:
                            print('Error - curve_fit failed for peak at {:}'.format(
                                    lines_in_region[i]))
                            print('Peak {:} added to output with intensity 0'.format(
                                    lines_in_region[i]))
                            
                            #Add in a 0 intensity entry for the failed peak
                            #Added on 1-7-19
                            temp = np.empty((1,4))
                            temp[0,0] = lines_in_region[i]
                            temp[0,1] = 0
                            temp[0,2] = 0
                            temp[0,3] = frame 
                            intens_arr = np.append(intens_arr,temp,axis=0)
                        except IndexError:
                            print('Error - No Intersections (peak {:} DNE)'.format(
                                    round(lines_in_region[i],5)))   
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
    #Added on 01-06-20 to incorporate spectrometer response calibration
    if ('correct' == True):
        for p in range(0,len(intensity_out[:,0])):
            #find the correction factor needed to apply the calibration
            #search in the calib (column 0) for nearest index of the entry which has the same
            #wavelength as the intensity_out array value
            calib_val = calib[find_nearest_index(calib[:,0],intensity_out[p,0]),1]
            intensity_out[p,1] = intensity_out[p,1] * 100 / calib_val
    if ((normalize == 'True') and (method == 'absolute')):
        max_val = np.amax(intensity_out[:,1])
        intensity_out[:,1] = scale_max * intensity_out[:,1] / max_val  
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
                print('No Peaks in wavelength range {:} : {:} - {:} nm'.format(countr,
                      regions[r,1],regions[r,2]))
    intensity_out = intensity_out_temp
    return intensity_out
