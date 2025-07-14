
import numpy as np
from synthetic_population import synthetic_population, total_number, number_density, ISO_diameters
from utils import true2ecc, ecc2mean, mean2tp, max_hc_distance_asteroid, year2sec, orb2cart, gal2ecl_cart, cart2orb, kepler, y2s
import matplotlib.pyplot as plt

from astropy.constants import au, GM_sun
population = 0

import time

## Constants
#au = 149597870700.0
#mu = 1.32712440042e20  # standard gravitational parameter of the Sun

time_of_simulation = 10.  # years
n0 = 1e-1  # number-density (for D > D_ref)
rm = 50

v_min = 100 # m/s
v_max = 1e5 # m/s
d_reff = 100


u_Sun = 1e4
v_Sun = 1.1e4
w_Sun = 7e3

sigma = [[1.2e4, 3.1e4, 2.6e4], [1.1e4, 2.3e4, 1.8e4],
         [0.9e4, 1.6e4, 1.5e4]]  # velocity dispersions for 3 stellar populations

vertex_deviation = [np.deg2rad(36), np.deg2rad(7), np.deg2rad(12)]  # vertex deviation for 3 stellar populations

"""
delimo populaciju po velicini i brzini zato sto od ova dva parametra zavisi odakle mogu da stignu do 
opservabilne sfere
"""

step = 100 # days


sigma_vx = sigma[0][population]
sigma_vy = sigma[1][population]
sigma_vz = sigma[2][population]
vd = vertex_deviation[population]




                
q, ecc, f, inc, Omega, omega = synthetic_population(rm = rm, n0 = n0,
                                                     v_min = v_min,
                                                     v_max = v_max,
                                                     u_Sun=u_Sun, v_Sun=v_Sun, w_Sun=w_Sun,
                                                     sigma_vx=sigma_vx, sigma_vy=sigma_vy,
                                                     sigma_vz=sigma_vz,
                                                     vd = vd, va = 0, R_reff = 696340000.,
                                                     speed_resolution=10, 
                                                     angle_resolution=60, dr=0.1)



a = q / (1 - ecc) * au.value
mm = np.sqrt(-GM_sun.value / a**3) 

M = np.zeros_like(q)
E = np.zeros_like(q)
for i in range(len(q)):
    E[i] = true2ecc(f[i], ecc[i])
    M[i] = ecc2mean(E[i], ecc[i])
    
r_init = a * (1 - ecc * np.cosh(E)) / au.value
    


    
'''
Propagating the population
'''
time_of_simulation = y2s(time_of_simulation)
step  =  step * 86400.

num = 0
time_propagation = 0.
current_number = []

while time_propagation < time_of_simulation:
    
    
    
    print(time_propagation / 86400)
    M_t = M + mm * time_propagation
    E_t = np.zeros_like(q)
    for i in range(len(q)):
        E_t[i] = kepler(ecc[i], M_t[i], accuracy = 1e-6)
    r = a * (1 - ecc * np.cosh(E_t)) / au.value
    current_number.append(np.sum(r < rm))

    num += 1
    time_propagation += step

  
  
flux = np.mean(np.diff(current_number) / (step))



n = 1e-4

n_total = n - flux * time_of_simulation / len(q) * n 


begining = time.time()
q, ecc, f, inc, Omega, omega, D = synthetic_population(rm = rm, n0 = n_total,
                                                     v_min = v_min,
                                                     v_max = v_max,
                                                     u_Sun=u_Sun, v_Sun=v_Sun, w_Sun=w_Sun,
                                                     sigma_vx=sigma_vx, sigma_vy=sigma_vy,
                                                     sigma_vz=sigma_vz,
                                                     vd = vd, va = 0, R_reff = 696340000.,
                                                     speed_resolution=10, 
                                                     angle_resolution=60, dr=0.1, d=[100, 60000], alpha=[-2.5])


print(f'duration = {time.time() - begining}')
a = q / (1 - ecc) * au.value
mm = np.sqrt(-GM_sun.value / a**3) 



initially_inside = int(n/n_total * len(q)) # which should be inside at the begining



M = np.zeros_like(q)
E = np.zeros_like(q)
for i in range(len(q)):
    E[i] = true2ecc(f[i], ecc[i])
    M[i] = ecc2mean(E[i], ecc[i])
    
# selektovanje onih koji treba da su unutra
    

    
    
M_lim = np.zeros_like(q)    
E_lim = np.arccosh(1/ecc * (1 - rm / (a / au.value)))

for i in range(len(q)):
    M_lim[i] = -ecc2mean(E_lim[i], ecc[i])
    
M_min = M_lim - mm * time_of_simulation

M_rand = (np.random.random(len(M_lim)) * (M_min - M_lim) + M_lim)
M_rand[:initially_inside] = M[:initially_inside]

E_rand = np.zeros_like(M_rand)



for i in range(len(E_rand)):
    E_rand[i] = kepler(ecc[i], M_rand[i], accuracy = 1e-6)
    


    
#r_init1 = a * (1 - ecc * np.cosh(E_rand)) / au.value






'''
Propagating the population (ulaze)
'''

num = 0
time_propagation = 0.
current_number1 = []

while time_propagation < time_of_simulation:
    
    
    
    print(time_propagation / 86400)
    M_t = M_rand[initially_inside:] + mm[initially_inside:] * time_propagation
    E_t = np.zeros_like(q[initially_inside:])
    for i in range(len(q[initially_inside:])):
        E_t[i] = kepler(ecc[initially_inside:][i], M_t[i], accuracy = 1e-6)
    r = a[initially_inside:] * (1 - ecc[initially_inside:] * np.cosh(E_t)) / au.value
    
    current_number1.append(np.sum(r < rm))

    time_propagation += step
    
    
    
'''
Propagating the population (ulaze)
'''

num = 0
time_propagation = 0.
current_number_in = []

while time_propagation < time_of_simulation:
    
    
    
    print(time_propagation / 86400)
    M_t = M_rand[:initially_inside] + mm[:initially_inside] * time_propagation
    E_t = np.zeros_like(q[:initially_inside])
    for i in range(len(q[:initially_inside])):
        E_t[i] = kepler(ecc[:initially_inside][i], M_t[i], accuracy = 1e-6)
    r = a[:initially_inside] * (1 - ecc[:initially_inside] * np.cosh(E_t)) / au.value
    
    current_number_in.append(np.sum(r < rm))
    
    
    
    
    num += 1
    time_propagation += step


plt.plot((np.array(current_number1) + np.array(current_number_in)) / initially_inside * 100)


    
