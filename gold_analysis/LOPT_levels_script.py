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
os.chdir(LOPT_folder)
output = ''
levels_dat = glob.glob('au_I_levels*.txt')
data_in = np.genfromtxt(levels_dat[0],delimiter='\t', dtype = str, skip_header = 1)
os.chdir(cwd)
level_labels = np.genfromtxt('table_auI_level_info.txt',skip_header = 1, delimiter ='\t',\ 
                             dtype = str)
for i in range(0,len(data_in[:,0])):    
    data_in[i,5] = (data_in[i,5]).split(",")[0]
#%% 
for i in range(0,len(level_labels[:,0])):
    output_temp = ''
    config = str(level_labels[i,1])
    total_j = level_labels[i,2]
    energy_l = str(round(float(level_labels[i,3]),2))
    energy_e = ''
    diff = ''
    uncert = ''
    no_lines = '0'
    if (i==0):
        energy_e = '0.0'
        diff = '0.0'
        uncert = '0.0'
        no_lines = data_in[0,5]
    if (i > 0):
        for j in range(1,len(data_in[:,0])):
            if (str(round(float(level_labels[i,0]))) == str(round(float(data_in[j,0])))):
                energy_e = str(round(float(data_in[j,1])))+str('~$\pm$~')+str(data_in[j,3])
                diff = round(float(data_in[j,1]) - float(energy_l),1)
                #uncert = data_in[j,3]
                no_lines = data_in[j,5]
    ref = str(level_labels[i,4])

    output_temp = str(config)+str('$_{')+str(total_j)+str('}$')+str('&')+str(energy_e)+str('&') \
    +str(energy_l)+str('&')+str(diff)+str('&') \ 
    + str(no_lines)+str('&')+str(ref)+str('\\\\')+str('\n')
    output = output + output_temp
            
text_file = open('AuI_levels_table.txt', "w")
text_file.write('%s' % output)
text_file.close()

os.chdir(paper_path)
text_file = open('AuI_levels_table.txt', "w")
text_file.write('%s' % output)
text_file.close()
os.chdir(cwd)            
    
    
            