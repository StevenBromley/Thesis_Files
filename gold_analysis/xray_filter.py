"""
Steve Bromley
10-30-19
X-Ray filter for CCD Data
"""
import os
import numpy as np
import csv
import glob
from numpy import genfromtxt
#Function for reading standard .csv form of .spe data
#I.e. stacked frames in 2 column format with
#Col 1 (Wavelength \ pixel) | Col 2 (Counts)
def read_prin(fil):
    f = open(fil)
    data = np.loadtxt(f,delimiter=',')
    inten = np.reshape(data[:,1],(-1,1340))
    return np.vstack([data[0:1340,0],inten])
############################################

#Max/min lambda for binning later; bounds will be adjusted depending on input
global_max_lambda = 0
global_min_lambda = 1000
#Threshold for an x-ray being real, i.e. an x-ray contaminated pixel
#has a ratio exceeding the threshold when compared to neighboring pixels
xray_thresh = 2.0
#Save main working directory:
cwd = os.getcwd()
#Number of bins for later:
data_path = 'C:/Users/steve_000/raw_data'
os.chdir(data_path)
manual_input = glob.glob("*.csv")
#Loop over ALL data files of interest in "manual_input" variable:
for data_in in manual_input: #glob.glob('*.csv'):
    os.chdir(data_path)   
    pin =  read_prin(data_in)
    #Total Number of Frames in data file:
    frames = len(pin[:,0])
    #Final peak output array from gaussian fits
    #Select range of frames for analysis:
    frame1 = 65
    frame2 = 110
    for frame in range(2,len(pin[:,0])-2):#range(frame1,frame2):
    ####################################
    #########      BEGIN       #########
    #########   X-RAY FILTER   #########
    #########  ON EACH FRAME   #########
    ####################################
        for k in range(0,3):
            for i in range(2,len(pin[0,:]) - 3):
#Check for each type of x-ray, i.e. 1 pixel, 2 pixel, 3 pixel
                inp = pin[frame,i]
                #2 pixel x-ray, or two side by side (common in later frames):
                if(((inp / pin[frame,i-1])>xray_thresh) and \
                   ((pin[frame,i+1]/pin[frame,i+2])>xray_thresh)):
                    pin[frame,i] = (pin[frame,i-1] + pin[frame,i+2]) / 2
                    pin[frame,i+1] = (pin[frame,i-1] + pin[frame,i+2]) / 2
                #3 pixel x-ray or 3 side by side (common toward end of shot):
                if ((inp/pin[frame,i-1] > xray_thresh) and \
                    (pin[frame,i+1]/pin[frame,i-1] > xray_thresh) and \
                    (pin[frame,i+2]/pin[frame,i+3] > xray_thresh)) :
                    pin[frame,i] = (pin[frame,i-1] + pin[frame,i+3]) / 2
                    pin[frame,i+1] = (pin[frame,i-1] + pin[frame,i+3]) / 2                                                            
                    pin[frame,i+2] = (pin[frame,i-1] + pin[frame,i+3]) / 2
                #1 pixel x-ray:
                #Check to see if x-ray is transient, i.e. an x-ray and not a peak:
                if((inp > pin[frame-1,i]) or (inp > pin[frame+1,i])):
                    #Remove x-ray by comparing to pixels on either side:
                    if(pin[frame,i]/pin[frame,i+1] > xray_thresh):
                        pin[frame,i] = (pin[frame,i+1] + pin[frame,i-1]) / 2
                    if(pin[frame,i]/pin[frame,i-1] > xray_thresh):
                        pin[frame,i] = (pin[frame,i+1] + pin[frame,i-1]) / 2
                    if (((inp/pin[frame-1,i]) > xray_thresh) and \
                        ((inp/pin[frame+1,i]) > xray_thresh)):
                        pin[frame,i] = (pin[frame-1,i] + pin[frame+1,i]) / 2
                    #check small, transient in time x-rays:
                    #Check based on change in amplitude in time:
                    if ((inp/pin[frame-1,i] > 1.4) and (inp/pin[frame+1,i] > 1.4)):
                        pin[frame,i] = (pin[frame-1,i] + pin[frame+1,i]) / 2
                    #We can compare to pixels on either side for smaller x-rays
                    if ((inp/pin[frame,i-1] > 1.4) and (inp/pin[frame,i+1] > 1.4)):
                        pin[frame,i] = (pin[frame,i-1] + pin[frame,i+1]) / 2        
    #Save as the same file name in the working directory:
    os.chdir(cwd)
    shot_num = re.findall("(\d+)",data_in)
    shot_num_int = int(shot_num[0])
    np.savetxt('{:}.csv'.format(shot_num_int), pin,delimiter =',')

