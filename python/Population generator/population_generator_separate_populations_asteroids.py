import os
import sys
import numpy as np
from synthetic_population import synthetic_population, total_number, number_density, ISO_diameters
from utils import true2ecc, ecc2mean, mean2tp, max_hc_distance_asteroid, year2sec, orb2cart, gal2ecl_cart, cart2orb, absolute_magnitude_asteroid

import time

number_of_jobs = int(sys.argv[1])
job_number = int(sys.argv[2])
population = int(sys.argv[3])



#number_of_jobs = 8
#job_number = 1


maximum_array_size=int(1e7/number_of_jobs)

# Constants
au = 149597870700.0
mu = 1.32712440042e20  # standard gravitational parameter of the Sun
epoch = 60796.00143922635 # epoch of the orbital elements (initial epoch of the survey simulation)

# Input parameters
maximum_array_size = int(1e5) # maximum number of objects per run
time_of_simulation = 10.  # years
n0 = 1e-1  # number-density (for D > D_ref)
d_reff = 100  # meters
v_min = 100 # m/s
v_max = 1e5 # m/s
d_min=10 # m
d_max = 100000 # m
alpha = -4 # maximum slope

V_cut = 24.5
albedo_max = 1.




division = 100 # speed range division


#resolution=[30, 60, 5.] # speed, angle, heliocentric distance
objects_per_file = int(1e5)



file_progress = 'A_' + stars + '_progress_job_'+ str(job_number) + '.txt'


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

# --------------------------------------------------------------

v_divided = np.linspace(v_min, v_max, division+1)


sigma_vx = sigma[0][population]
sigma_vy = sigma[1][population]
sigma_vz = sigma[2][population]
vd = vertex_deviation[population]


    
colors = [[2.31, 0.73, -0.33, -0.50, -0.68, 0.15]]
sfd = np.arange(1.0, 4.2, 0.2)
albedo = np.array([1.0])
N_max=n0/number_of_jobs*((d_reff/d_min)**(-alpha) - (d_reff/d_max)**(-alpha)) 
start_numbers = np.zeros([len(sfd), len(albedo)], dtype = 'int')

for slope in sfd:     
    
    orbit_file = 'results_new/asteroids/A_' + str(slope) + '_orbit_file_' + str(job_number) + '.txt'
    pp_file = 'results_new/asteroids/A_' + str(slope) + '_pp_file_' + str(job_number) + '.txt'

    if os.path.exists(orbit_file):
        os.remove(orbit_file)
    if os.path.exists(pp_file):
        os.remove(pp_file)

    file_writer_orbit = open(orbit_file, 'ab')
    np.savetxt(file_writer_orbit, ['ObjID,FORMAT,q,e,inc,node,argPeri,t_p_MJD_TDB,epochMJD_TDB'], fmt='%s')

    file_writer_pp = open(pp_file, 'ab')
    np.savetxt(file_writer_pp, ['ObjID,H_g,u-r,g-r,i-r,z-r,y-r,GS'], fmt='%s')
          
    file_writer_orbit.close()
    file_writer_pp.close()

# =============================================================================
#                               SIMULATION
# =============================================================================


"""
Maximum heliocentric distance where the largest object from the population can be observed
# The maximum diameter is set to 1 km
"""

hc_max = max_hc_distance_asteroid(d_max, albedo_max, V_cut)  # This is OK (double checked!)

"""
Since the OIF simulation lasts some predefined time (10 year in our case), objects which are initially
further away from the Sun than hc_max may reach this heliocentric distance during the simulation time. 
This means we need to increase the model sphere a bit. This increment is how much the fastest object 
from the population can travel during that time, or 1 year * v_max 

This gives the radius of our model sphere rm = hc_max + 1 year * v_max.

Object inside this sphere might be at heliocentric distance where they are observable or
they might be further away but can reach the observable heliocentric distance during the simulation time.
Objects outside this are for sure not observable during the simulation time. 

If we increase simulation time (e.g. to 10 years), then object which are initially very far away from the Sun
can maybe reach observable heliocentric distance during that time. This increases the model sphere, number of 
objects, computational resources....
"""

"""
Output files. In order to avoid very big files, we write maximum 20M object in one file (~1 GB)
"""

file_number=0

br=0 # total number of final objects



"""
sada za svako v racunamo rm
"""

rm = np.zeros(division)
for k in range(division):
    rm[k] = hc_max + year2sec(time_of_simulation) * v_divided[k+1] / au  # This is OK (double checked!)

# total number-density which iuncludes all objects within the defined size range   
n_total=number_density(n0/number_of_jobs, d_ref=d_reff, d=[d_min, d_max], alpha=[alpha])[-1]
    
pocetak=time.time()

iso_total=np.zeros_like(rm)



for k1 in range(len(rm)):
    
    pocetak1=time.time()
    
    ISO_total_number = total_number(rm=rm[k1], n0=n_total, v_min=v_divided[k1], v_max=v_divided[k1 + 1],
                                    u_Sun=u_Sun, v_Sun=v_Sun, w_Sun=w_Sun,
                                    sigma_vx=sigma_vx, sigma_vy=sigma_vy, sigma_vz=sigma_vz,
                                    vd = vd, va = 0, R_reff = 696340000.,
                                    speed_resolution=10, 
                                    angle_resolution=60, dr=0.5,
                                    d_ref=d_reff, maximum_array_size=maximum_array_size)  # This is OK (double checked)

           
    if ISO_total_number>0:
        number_of_runs = 1

        if ISO_total_number > maximum_array_size:
            number_of_runs = int(np.ceil(ISO_total_number / maximum_array_size))
            
        run=0

        while run<number_of_runs:
            
            print(run)
            
            try:
                q_out = np.array([])
                e_out = np.array([])
                f_out = np.array([])
                inc_out = np.array([])
                Omega_out = np.array([])
                omega_out = np.array([])
    
                np.savetxt(file_progress,
                           ['Run number {} out of {}, v={}.'.format(
                               run + 1, number_of_runs, [v_divided[k1], v_divided[k1 + 1]])], fmt='%s')
    
                print(
                    '\n ------------------------------ \n Run number {} out of {}, v={}.'.format(
                        run + 1, number_of_runs, [v_divided[k1], v_divided[k1 + 1]]))
    
                print('total number=', ISO_total_number)
                
                

                
                q, ecc, f, inc, Omega, omega = synthetic_population(rm=rm[k1], n0=n_total/number_of_runs,
                                                                     v_min=v_divided[k1],
                                                                     v_max=v_divided[k1 + 1],
                                                                     u_Sun=u_Sun, v_Sun=v_Sun, w_Sun=w_Sun,
                                                                     sigma_vx=sigma_vx, sigma_vy=sigma_vy,
                                                                     sigma_vz=sigma_vz,
                                                                     vd = vd, va = 0, R_reff = 696340000.,
                                                                     speed_resolution=10, 
                                                                     angle_resolution=60, dr=0.5,
                                                                     d_ref=d_reff, maximum_array_size=maximum_array_size)  # This is OK (double checked)
                
                if len(q) > 0:
                
                    for enum_slope, slope in enumerate(sfd):
                        
                        N_ref = n0/number_of_jobs*((d_reff/d_min)**slope - (d_reff/d_max)**slope) 
                        number = int(N_ref/N_max * len(q))
            
                        selection = np.random.randint(0, len(q), number)
                        
                        q1 = q[selection]
                        ecc1 = ecc[selection]
                        inc1 = inc[selection]
                        Omega1 = Omega[selection]
                        omega1 = omega[selection]
                        f1 = f[selection]
                        n, nt = number_density(n0/number_of_jobs, d_reff, [d_min, d_max], [-slope])
                        D = ISO_diameters(n, [d_min, d_max], [-slope], len(q1))
                    
                        hc_max = max_hc_distance_asteroid(D, albedo, 24.5)
                        
                        selection1 = q1 < hc_max  # This is OK (double checked)
            
                        ecc2 = ecc1[selection1]
                        f2 = f1[selection1]
                        inc2 = inc1[selection1]
                        Omega2 = Omega1[selection1]
                        omega2 = omega1[selection1]
                        q2 = q1[selection1]
                        D2 = D[selection1]
                        hc_max2 = hc_max[selection1]
             
                        '''
                        SELECTION No. 2
                        Now we check if the object is inside observable time at initial moment
                        Equation of hyperbolic orbit:
                        r = a*(1-e*cosh(E)), where E is hyperbolic anomaly
                        
                        From this equation, given eccentricity, semi-major axis and maximum observable heliocentric distance
                        we can calculate critical hyperbolic anomaly (when ISO is exactly at hc_max) 
                        '''
                        Ecr = np.arccosh(1 / ecc2 - hc_max2 / ecc2 / (q2 / (1 - ecc2)))  # This is OK (double checked) this is always positive
        
                        """
                        corresponding critical mean anomaly (maximum where an object of a given size can be observed)
                        from hyperbolic Kepler equation
                        M = e *sinh(E) - E
                        """
        
                        M_max = ecc2 * np.sinh(Ecr) - Ecr  # This is OK (double checked) this is always positive
        
                        """
                        We calculate M for every object from the population
                        """
                        
                        mean_motion = np.sqrt(mu / (np.abs(q2 / (1 - ecc2)) * au) ** 3)
                        
                        M = np.zeros(len(q2))  # current mean anomaly
                        for i in range(len(q2)):
                            E = true2ecc(f2[i], ecc2[i])  # This is OK (double checked)
                            M[i] = ecc2mean(E, ecc2[i])  # This is OK (double checked)
        
                        """
                        We calculate mean motion
                        """
                          # mean motion (This is OK) (double checked)
        
                        """
                        Finally, we calculate minimum mean anomaly from which an object can reach hc_max during the simulation time
                        """
                        M_min = -M_max - mean_motion * year2sec(time_of_simulation)  # This is OK (double checked)
        
                        """                        
                        If object is outside observable sphere and it is on outgoing branch or orbit (Mean anomaly > M_max) it is excluded
                        because it is surely non observable. For objects with positive mean anomaly (those which are on the outgoing branch)
                        we take only those which are inside their observable spheres
                        
                        
                        For objects on incoming branches we take those which are currently observable but also those whose mean anomaly 
                        is larger than M_min. This means that those object will reach their observable sphere during the OIF simualtion time. 
                        """
        
                        selection2 = np.logical_and(M > M_min, M < M_max)  # This is OK (double checked)
                        e_out = ecc2[selection2]
                        f_out = f2[selection2]
                        inc_out = inc2[selection2]
                        Omega_out = Omega2[selection2]
                        omega_out = omega2[selection2]
                        q_out = q2[selection2]
                        D_out = D2[selection2]
                        H_out = absolute_magnitude_asteroid(D_out, a)
                    
                        for i in range(len(q_out)):
                            x,y,z,vx,vy,vz=orb2cart(omega_out[i], Omega_out[i], inc_out[i], e_out[i], q_out[i]/(1-e_out[i])*au, 0.1, mu) # we put E=0.1 (just not at the pericenter) because it does not impact the conversion of o, O, inc (This is OK - double checked)
                            
                            x,y,z=gal2ecl_cart(x,y,z) # This is OK (double checked)
                            vx,vy,vz=gal2ecl_cart(vx,vy,vz)  # This is OK (double checked)  
                            omega_out[i], Omega_out[i], inc_out[i]=cart2orb(x,y,z,vx,vy,vz, mu)[:3] # This is OK (double checked)
                            
                            
                        tp_out=np.zeros(len(q_out))
    
                        for i in range(len(q_out)):
    
                            tp_out[i] = mean2tp(ecc2mean(true2ecc(f_out[i], e_out[i]), e_out[i]), q_out[i] / (1 - e_out[i]), epoch)  # This is OK (double checked)
                           
                        start_number = start_numbers[enum_slope][enum_a]
                        end_number = start_number + len(q_out)
                        
                        prefix = 'A_'  + stars + '_' + str(slope) + '_' + str(a) + '_'
                        
                        names = (np.arange(start_number, end_number, dtype = 'int')).astype(str)
                        
                        full_names = [prefix + x for x in names]
                        

                        full_names = np.array(full_names, dtype='str')
                        
                        COM = np.array([['COM']] * len(q_out), dtype='str')
        
                        orbit_file = 'results_new/asteroids/A_'  + stars + '_' + str(slope) + '_' + str(a) + '_orbit_file_' + str(job_number) + '_.txt'
                        pp_file = 'results_new/asteroids/A_'  + stars + '_'  + str(slope) + '_' + str(a) + '_pp_file_' + str(job_number) + '_.txt'
                        file_writer_orbit = open(orbit_file, 'ab')
                        file_writer_pp = open(pp_file, 'ab')
                        
        
                        # writing file for Sorcha simulation
                        np.savetxt(file_writer_orbit, np.column_stack([full_names, COM, np.round(q_out, 6), np.round(e_out, 6), np.round(inc_out, 3),
                                 np.round(Omega_out, 3), np.round(omega_out, 3), np.round(tp_out, 6), np.round(np.ones(len(q_out)) * epoch, 6)]),
                                       fmt='%s', delimiter = ',')
            
            
                        np.savetxt(file_writer_pp, np.column_stack([full_names, np.round(H_out, 2), colors * len(q_out)]),
                                       fmt='%s', delimiter = ',')
            
                        start_numbers[enum_slope] = end_number
                    
                    
                        file_writer_orbit.close()
                        file_writer_pp.close()
        
        
        
            except Exception as e: # if anything goes wrong, just repeat as if nothing happened
                print(f"An error occurred: {e} - asteroids")
             
                run=run-1
                
            run=run+1
              


trajanje=time.time()-pocetak

np.savetxt(file_progress, [trajanje], fmt='%s')
file_progress.close()
