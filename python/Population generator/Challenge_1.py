
import numpy as np
from synthetic_population import synthetic_population, total_number, number_density, ISO_diameters
from utils import true2ecc, ecc2mean, mean2tp, max_hc_distance_asteroid, year2sec, orb2cart, gal2ecl_cart, cart2orb, kepler, y2s
import matplotlib.pyplot as plt

from astropy.constants import au, GM_sun
population = 0

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
time = 0.
current_number = []

while time < time_of_simulation:
    
    print(time / 86400)
    M_t = M + mm * time
    E_t = np.zeros_like(q)
    for i in range(len(q)):
        E_t[i] = kepler(ecc[i], M_t[i], accuracy = 1e-6)
    r = a * (1 - ecc * np.cosh(E_t)) / au.value
    
    
    
    current_number.append(np.sum(r < rm))
    
    num += 1
    time += step
#    break
    
    
flux = np.mean(np.diff(current_number) / (step) )

n = 1e-4

n_total = n - flux * time_of_simulation / len(q) * n 

q1, ecc1, f1, inc1, Omega1, omega1, D1 = synthetic_population(rm = rm, n0 = n_total,
                                                     v_min = v_min,
                                                     v_max = v_max,
                                                     u_Sun=u_Sun, v_Sun=v_Sun, w_Sun=w_Sun,
                                                     sigma_vx=sigma_vx, sigma_vy=sigma_vy,
                                                     sigma_vz=sigma_vz,
                                                     vd = vd, va = 0, R_reff = 696340000.,
                                                     speed_resolution=10, 
                                                     angle_resolution=60, dr=0.1, d=[100, 60000], alpha=[-2.5])


a1 = q1 / (1 - ecc1) * au.value
mm1 = np.sqrt(-GM_sun.value / a1**3) 



fraction = int(n/n_total * len(q1))









M1 = np.zeros_like(q1)
E1 = np.zeros_like(q1)
for i in range(len(q1)):
    E1[i] = true2ecc(f1[i], ecc1[i])
    M1[i] = ecc2mean(E1[i], ecc1[i])
    
# selektovanje onih koji treba da su unutra
    

    
    
M_lim = np.zeros_like(q1)    
E_lim = np.arccosh(1/ecc1 * (1 - rm / (a1 / au.value)))

for i in range(len(q1)):
    M_lim[i] = -ecc2mean(E_lim[i], ecc1[i])
    
M_min = M_lim - mm1 * time_of_simulation

M_rand = (np.random.random(len(M_lim)) * (M_min - M_lim) + M_lim) * 0.96
M_rand[:fraction] = M1[:fraction]

E_rand = np.zeros_like(M_rand)



for i in range(len(E_rand)):
    E_rand[i] = kepler(ecc1[i], M_rand[i], accuracy = 1e-6)
    


    
r_init1 = a1 * (1 - ecc1 * np.cosh(E_rand)) / au.value






'''
Propagating the population (ulaze)
'''

num = 0
time = 0.
current_number1 = []

while time < time_of_simulation:
    
    print(time / 86400)
    M_t1 = M_rand[fraction:] + mm1[fraction:] * time
    E_t1 = np.zeros_like(q1[fraction:])
    for i in range(len(q1[fraction:])):
        E_t1[i] = kepler(ecc1[fraction:][i], M_t1[i], accuracy = 1e-6)
    r1 = a1[fraction:] * (1 - ecc1[fraction:] * np.cosh(E_t1)) / au.value
    
    
    
    current_number1.append(np.sum(r1 < rm))
    
    num += 1
    time += step
    
    
    
'''
Propagating the population (ulaze)
'''

num = 0
time = 0.
current_number_in = []

while time < time_of_simulation:
    
    print(time / 86400)
    M_t1 = M_rand[:fraction] + mm1[:fraction] * time
    E_t1 = np.zeros_like(q1[:fraction])
    for i in range(len(q1[:fraction])):
        E_t1[i] = kepler(ecc1[:fraction][i], M_t1[i], accuracy = 1e-6)
    r1 = a1[:fraction] * (1 - ecc1[:fraction] * np.cosh(E_t1)) / au.value
    
    
    
    current_number_in.append(np.sum(r1 < rm))
    
    num += 1
    time += step


plt.plot((np.array(current_number1) + np.array(current_number_in)) / len(r_init) * 100)


    
