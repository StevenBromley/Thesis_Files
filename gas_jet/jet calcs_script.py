# -*- coding: utf-8 -*-
"""
SJB
09-30-19
Gas Jet Code
"""

import numpy as np
import math
from scipy import integrate
import scipy
import matplotlib.pyplot as plt

def find_nearest(array, value):
    idx = (np.abs(array - value)).argmin()
    return idx

#%%
   
data_in = np.genfromtxt('jet_inp.csv',delimiter = ',', skip_header = 1)
n_out = np.zeros((4,5))
#n_out format:
#Col 1   Col 2   Col 3   Col 4   Col 5
#N       n_o      Q1      Q2     Qtot

#%%
k_b = 1.38064852e-23 #m^2 kg s^-2 K^-1
#M = 6.6464764e-27#He

#N2
M = 4.651737683e-26 #kg

#Nozzle geometry:
r_noz = .007 * .0254 # (radius (in) * conversion to meter)
A_nozz = np.pi * (r_noz)**2
#current:
sep = .57 * .0254 # .57 in times conversion to meters
#sep = 0.25 * 0.0254
#new:
#sep = .25 * .0254
z_skim = sep

#new design: 
#sep = .2 *.0254

#Inputs:
T_o = 300#kelvin
T_T = 0.75 * T_o

P_B = 9e-3
#Old design: 1mbar (1oo pA). If we use the hydrgen pump can we do 10mbar?)
#P_o = 50000   #Pascals. 100 Pa = 1 mbar; new design at 100 mbar
P_o = 500
n_o = P_o / (k_b * T_o)
#current design:
delta_r = .012 * 0.0254 #in to meter)

#Mach Disk location:
mach = 2 * r_noz * 0.67 * math.sqrt(P_o / P_B)

#new design:
#delta_r = 0.05 * .0254
#Number leaving the nozzle per second:
N_nozz = (3*P_o*A_nozz / 16) * math.sqrt(15/(k_b*T_o*M))
#Jet Velocity
v_jet = math.sqrt(5*k_b * T_o / M)
#%Maximum allowed transverse temp:
v_perp = (v_jet * (delta_r)/ z_skim)
T_max = M * (v_perp * v_perp) /(k_b*T_o)



#v distribution; x = velocity here

#3d Maxwellian @ T_nozz (Tzz = 3/4 To)
v_func = lambda x: 4*np.pi * x**2 * (M/(2*np.pi*k_b*T_T))**(3/2) * math.e**(-M*x**2 / (2*k_b*T_T))

#integrate # that get past skimmer
#chi = scipy.integrate.quad(v_func,0,v_perp)[0]
chi = scipy.integrate.quad(v_func,0,v_perp)[0]
print('Percent passing through skimmer: {:}'.format(100*chi))
#%%
#Divergence Half angle
alpha = np.arctan(1/math.sqrt(chi))
print(alpha)
#%%
#on-axis density as a function of distance
z = np.linspace(0.0001,1,1000)

#%%
N_skim = delta_r * delta_r * chi * N_nozz/ (sep * sep)

A_z = np.pi * z**2 * (1/chi)
z = 16 * .0254
#z = .2
nz = (n_o * 3 * math.sqrt(3)/16 * A_nozz * chi/ (np.pi * (z*z))) * 1e-6 #cm^-3

#Currently:
#z_int = 0.188 #meters
#z_int = .075
z_int = 0.4

z_skim = .57*.0254 #converted to meters
nz_int = (n_o * 3 * math.sqrt(3)/16 * A_nozz * chi/ (np.pi * (z_int*z_int)))
nz_int_skim = (n_o * 3 * math.sqrt(3)/16 * A_nozz * chi/ (np.pi * (z_skim*z_skim)))

nz_skim_arr = np.empty((0,2))
temp = np.empty((1,2))
for i in range(0,len(data_in[:,0])):
    N_skim_model = (3*data_in[i,0] * 133.322 *A_nozz / 16) * math.sqrt(15/(k_b*T_o*M)) * chi
    temp[0,0] = data_in[i,0]
    temp[0,1] = N_skim_model
    nz_skim_arr = np.append(nz_skim_arr,temp,axis=0)


#%%
#%%
####Calculate throughput
#Pump Speeds
Se1 = 33.11
Se2 = 6.01
T = T_T
M = 2.325867e-26 * 2
#nozzle diam =0.0135in, ra
#rad = (3.429e-4 / 2)
rad = 3.048e-4 #radius of skimmer (0.012 in.) in meters
A = np.pi * rad * rad #area of skimmer opening
#Inputs:
#Calculate density inside skimmer entrance from the pressure changes:
for i in range(1,len(data_in[:,0])):
    P_PKR = data_in[i,2] - data_in[0,2]
    deltaP = data_in[i,1] - data_in[i,2] * 0.133322368 #converts to Pa
    Q2 = Se2 * P_PKR #L/s times torr 
    N = Q2 * 0.133322368 / (k_b*T) # L/s times Torr x (torr to Pa) x (L to cm^3) / kT
    n_o =  N / (1e4 * A * v_jet * 1e2) #last term converts to cm^-3
    leak = 0.05 * deltaP / (k_b * T)
    n_out[i-1,0] = N
    n_out[i-1,2] = leak
    n_out[i-1,1] = n_o
    n_out[i-1,3] = Q2

#%%
p_lin = np.linspace(0,1000, 100000)
model_num = (3*p_lin[:] * 133.322 *A_nozz / 16) * math.sqrt(15/(k_b*T_o*M)) * chi

#%%
plt.clf()
plt.figure(figsize = (10,6))
plt.grid(True)
plt.xlim(0,10)
#plt.xlim(1.5,6)
plt.rcParams.update({'font.size': 14})
#plt.title('Number of particles per second exiting skimmer')
plt.yscale('log')
plt.ylabel('Particles leaving skimmer [s$^{-1}$]', fontsize = 16)
plt.xlabel('Nozzle pressure [Torr]', fontsize = 16)
#plt.plot(data_in[1:,0],n_out[:,0], marker = 'o', linestyle = '--', color = 'k', label = 'Total Gas Load')
plt.plot(data_in[1:,0],abs(n_out[:,0]), marker = 'x', linestyle = '--', color = 'k', label = 'Experimental Test') 
#plt.plot(data_in[1:,0],n_out[:,2], marker = 's', linestyle = '--', color = 'b', label = 'Estimated Skimmer Leakage')
plt.scatter(nz_skim_arr[1:,0], nz_skim_arr[1:,1], label = 'Model', marker = 'o', color = 'red')
plt.plot(p_lin,model_num, color = 'red')
plt.legend()

#%%
test = n_out[:,0] / nz_skim_arr[1:,1]

#%%
#print('Interaction Density: {:e}'.format(nz * 1e-6))
print('Skimmer Num. per sec: {:e}'.format(N_skim))
print('T perp: {:} K'.format(T_max))
print('Divergence: {:}'.format(alpha))
print('Interaction Region n(z): {:e}'.format (nz_int * 1e-6))
print('Skimmer Exit Density: {:e}'.format(nz_int_skim * 1e-6))
#y = z*z
#n_z = n_at_z(n_o,A_nozz,chi,z)
"""
plt.clf()
plt.plot(z[:],nz)
plt.yscale('log')
#plt.plot(z,A_z)
#plt.ylim(1e14,1e17)
plt.xlabel('z (m)')
plt.ylabel('Central Density (#/m^3)')
plt.show()
"""