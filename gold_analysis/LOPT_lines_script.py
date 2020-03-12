# -*- coding: utf-8 -*-
"""
SJB
06-17-19
script for turning data files into latex text
"""
#######################################
import sys
#define path for steve_lib.py
sys.path.insert(0, 'C:/Users/steve_000/Desktop/Research_Files/PRA_gold_spectroscopy')
#IMPORT EVERYTHING
from steve_lib import *
import math
from math import floor
def round_to_1(x):
    return round(x, -int(floor(math.log10(abs(x)))))
#######################################
paper_path = 'C:/Users/steve_000/Desktop/paper_path'
LOPT_folder = 'C:/Users/steve_000/Desktop/AuI_LOPT'
cwd = os.getcwd()

#%%
os.chdir(LOPT_folder)
au_I_lines_dat = glob.glob('au_I_lines*.txt')
data_in = np.genfromtxt(au_I_lines_dat[0],delimiter ='\t', dtype = str,skip_header=1)
levels_dat = glob.glob('au_I_levels*.txt')
levels = np.genfromtxt(levels_dat[0],delimiter='\t', dtype = str, skip_header = 1)
os.chdir(cwd)
level_labels = np.genfromtxt('table_auI_levels_info.txt',skip_header = 1, \ 
                             delimiter ='\t', dtype = str)
shot_list = genfromtxt('shots.csv',delimiter = ',')
references = np.genfromtxt('AuI_line_refs.txt', delimiter ='\t', usecols = (0,3),dtype =str)
#Reference File format:
#Wavelength (nm) | Lower level | Upper Level | Reference Text
#I use level labels according to the scheme in the level_labels table
#ex: 454 1 20 PS means the line wavelength is 454 nm from level 20 to level 1, and the displayed text
#in the compiled table will read 'PS'
#Grab the regions w/ wavelength ranges:
#Spectrometer wavelength regions
regions = genfromtxt('regions_with_wavelengths_AuI.txt', delimiter ='\t')
#Comments to put in the table, e.g. 'Blended Line' etc.
#Comments use 2-col format: Wavelength (nm) | Comment Text
comments = np.genfromtxt('AuI_line_comments.txt', delimiter ='\t', dtype = str, usecols = (0,1))

#%%Remove commas from the number of lines part:
for i in range(0,len(levels[:,0])):    
    levels[i,5] = (levels[i,5]).split(",")[0]
#%%grab relative intensities
line_list = np.empty((0,1))
for i in range(0,len(data_in[:,0])):
    if (data_in[i,1] != '_'):
        line = float(data_in[i,1])
        line = line / 10
        line_list = np.append(line,line_list)
line_list = line_list[line_list[:].argsort()]        
probe2 = 'C:/Users/steve_000/Desktop/xray_filtered_data'
#%%
output = ''
#%%
#Choose whether you want the Ritz wavelengths to come from the experimental OR literature energies
#ritz_source = 'theory'
ritz_source = 'experimental'

line_intensities = intensity_search_v2(line_list,frame_in=83,normalize = 'True',scale_max = 100,\ 
                                       method = 'relative', regions = regions,\ 
                                       data_directory=probe2,frame_source='file')
for r in range(0,len(regions[:,0])):  
    #line_intensities = np.zeros((2,10))
    for i in range(0,len(data_in[:,0])):
        if (data_in[i,1] != '_'):
            if (regions[r,1] < float(data_in[i,1]) / 10 < regions[r,2]):
                output_temp = ''
                stdev = str(round_to_1(float(data_in[i,2]) / 10))
                wavelength = str(round(float(data_in[i,1]) / 10,2))+str('(') \ 
                +stdev.split(sep='.')[1]+str(')')
                for j in range(0,len(level_labels[:,0])):
                    if (data_in[i,13] != 'G0'):
                        if (str(round(float(data_in[i,13]))) == level_labels[j,0]):
                            lower_level = str(level_labels[j,1])
                            lower_j = str(level_labels[j,2])
                            lower_e = float(level_labels[j,3])
                    else:
                        lower_level = str(level_labels[0,1])
                        lower_j = str(level_labels[0,2])
                        lower_e = float(0)
                        
                    if (str(round(float(data_in[i,14]))) == level_labels[j,0]):
                        upper_level = str(level_labels[j,1])
                        upper_j = str(level_labels[j,2])
                        upper_e = float(level_labels[j,3])
                
                if (ritz_source == 'theory'):
                    ritz = val_vac_to_air(1e8 / (upper_e - lower_e)) / 10
                    #ritz = round(0.9997e7 / (upper_e - lower_e),2)
            
                if (ritz_source == 'experimental'):
                    ritz = float(data_in[i,5]) / 10
            
                #keep sig figs for ritz:
                diff = str(round(float(data_in[i,1])/10 - float(ritz),2))
                
                #round ritz for table
                ritz = str(round(float(ritz),2))
                ###########################   INTENSITY   ########
                check = 0
                for k in range(0,len(line_intensities[:,0])):
                    if (abs((float(data_in[i,1]) / 10) - line_intensities[k,0]) < .01):
                        #intens = str('{:.2E}'.format(round(line_intensities[k,1],3)))
                        intens = str('{:}'.format(int(round(line_intensities[k,1]+.3,0))))
                        check = 1
                if (check == 0):
                        intens = str('1*')            
                ###########################   REFERENCE   ########
                ref = ''
                for k in range(0,len(references[:,0])):
                    if (abs((float(data_in[i,1]) / 10) - float(references[k,0])) < 0.02):
                        ref = references[k,1]
                if (ref == ''):
                    ref = 'New'
                comment = ''
                for p in range(0,len(comments[:,0])):
                    if (abs((float(data_in[i,1]) / 10) - float(comments[p,0])) < 0.03):
                        comment = comments[p,1]
                #################### OUTPUT ######################
                output_temp = str(intens) + '&'+str(wavelength) + str('&') + str(ritz) \ 
                + str('&') + str(round(float(diff),2)) + str('&') + str(lower_level) \ 
                + str('&') +  str(lower_j) + str('&') + str(upper_level) + str('&') \ 
                + str(upper_j) + str('&') + str(comment) + str('&') +str(ref) \ 
                + str('\\\\') + str('\n') 
                #Note that python interprets the \f string as form feed and
                #NOT a literal '\f' in text.
                output = output + output_temp
            
            #############   WAVELENGTH REGIONS  ##############
    if ((r != 16) or (r != 15)):
        output = output + str('\cline{1-10}') + str('\n')

text_file = open('AuI_lines_table.txt', "w")
text_file.write('%s' % output)
text_file.close()

os.chdir(paper_path)
text_file = open('AuI_lines_table.txt', "w")
text_file.write('%s' % output)
text_file.close()
os.chdir(cwd)              
            