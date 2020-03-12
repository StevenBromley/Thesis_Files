# -*- coding: utf-8 -*-
"""
Calculating Ritz Wavelengths using Maria (1997) level energies for Au II
and comparing to an observed line list
SJB
"""

import sys
#define path for steve_lib.py
sys.path.insert(0, 'C:/Users/steve_000/Desktop/Research_Files/PRA_gold_spectroscopy')
#IMPORT EVERYTHING
from steve_lib import *

def frac_to_float(arr):
    temp = np.empty((1,1))
    output = np.empty((0,1))
    for i in range(0,len(arr)):
        temp = float(Fraction(arr[i]))
        output = np.append(output,temp)
    return output
    
def wavenumber_to_nm(val):
    nm = 1e9 * 1.98644568e-25 / (val * 1.9863E-23)
    return(nm)


def convert_to_float(frac_str):
    try:
        return float(frac_str)
    except ValueError:
        num, denom = frac_str.split('/')
        try:
            leading, num = num.split(' ')
            whole = float(leading)
        except ValueError:
            whole = 0
        frac = float(num) / float(denom)
        return whole - frac if whole < 0 else whole + frac    
    
#%Load in the master level file:
all_levels = genfromtxt('table_auII_levels_info.txt', delimiter ='\t', skip_header = 1, usecols = (0,1,2,3), dtype = str)        
for i in range(0,len(all_levels[:,0])):
    all_levels[i,1] = all_levels[i,1].replace('"','') 

temp_odd = np.empty((1,3))
odd = np.empty((0,3))
temp_even = np.empty((1,3))
even = np.empty((0,3))
#sort the master level file. In my case I've denoted the odd levels by 'circ'
#in the label in the above 'all_levels' variable
#The file 'table_auII_levels_info' is used for other scripts also.
for i in range(0,len(all_levels[:,0])):
    if (all_levels[i,2] != ''):
        if ('circ' in all_levels[i,1]):
            temp_odd[0,0] = all_levels[i,0]
            temp_odd[0,1] = all_levels[i,2]
            temp_odd[0,2] = all_levels[i,3]
            odd = np.append(odd,temp_odd,axis=0)
        else:
            temp_even[0,0] = all_levels[i,0]
            temp_even[0,1] = all_levels[i,2]
            temp_even[0,2] = all_levels[i,3]
            even = np.append(even,temp_even,axis=0)

#Load in a list of lines to compare to. Ex: my 2nd_6cm gold line list
lines = genfromtxt('2nd_6cm_full.txt', delimiter ='\t')
#%%
temp = np.empty((1,8))
output = np.empty((0,8))
#Loop over all possible combinations of odd/even levels
for i in range(0,len(even[:,0])):
    for j in range(0,len(odd[:,0])):
        #Enforce Delta J = 0, +/- 1
        if (abs(even[i,1] - odd[j,1]) <= 1):
            #If Delta J is satisfied, calculate the Ritz wavelength.
            #Also note that the val_vac_to_air function input MUST be in angstroms, so I've done
            #that conversion here and converted back to nm after running the ritz wavelength through
            #the function
            ritz = val_vac_to_air(10 * wavenumber_to_nm(abs(even[i,2] - odd[j,2]))) / 10
            #loop over the lines file. If any lines are within the threshold (0.06 nm), save all the 
            #information about that line:
            for k in range(0,len(lines[:,0])):
                if (abs(lines[k,0] - ritz) < .06):
                    if (odd[j,2] > even[i,2]):
                        temp[0,0] = lines[k,0] 
                        temp[0,1] = ritz
                        temp[0,2] = even[i,2]
                        temp[0,3] = even[i,1]
                        temp[0,4] = odd[j,2]
                        temp[0,5] = odd[j,1]
                        temp[0,6] = even[i,0]
                        temp[0,7] = odd[j,0]
                        output = np.append(output,temp,axis=0)
                    else:
                        temp[0,0] = lines[k,0] 
                        temp[0,1] = ritz
                        temp[0,2] = odd[j,2]
                        temp[0,3] = odd[j,1]
                        temp[0,4] = even[i,2]
                        temp[0,5] = even[i,1]
                        temp[0,6] = odd[j,0]
                        temp[0,7] = even[i,0]
                        output = np.append(output,temp,axis=0)
                    
#Save
output = output[output[:,0].argsort()]        
np.savetxt('Au_II_matches.txt', output, delimiter ='\t')
                