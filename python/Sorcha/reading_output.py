import pandas as pd
import numpy as np
import matplotlib.pyplot as plt



slope_values = np.arange(1, 4.2, 0.2)

cadence_file = 'cadence_file/baseline_v4.3.5_10yrs.db'

inc_detected = []
ratio = []
mean_inc = []
broj = []
for slope in slope_values:
    slope_formated = f"{round(slope, 1):.1f}"
    try:
        inc_temp = []
        for i in range(1, 21):
            
            orbit_file = 'asteroids/A_M_' + slope_formated + '_1.0_orbit_file_' + str(i) + '_.txt'
            pp_file = 'asteroids/A_M_' + slope_formated + '_1.0_pp_file_' + str(i) + '_.txt'
            output_file = 'results/no_trailling/A_M_' + slope_formated + '_1.0_output_file_' + str(i) + '.csv'
            
            output = pd.read_csv(output_file)
            
            inc = np.loadtxt(orbit_file, skiprows = 1, delimiter = ',', usecols = (4), unpack = True)
            ID_all = np.loadtxt(orbit_file, skiprows = 1, delimiter = ',', usecols = (0), unpack = True, dtype = str)
            ID_detected = list(set(output['ObjID']))
            
            for j in range(len(ID_detected)):
                inc_temp.append(inc[ID_all == ID_detected[j]][0])
                
        mean_inc.append(np.mean(inc_temp))
        inc_detected.append(inc_temp)
        
        ratio.append(np.sum(np.array(inc_temp) > np.pi/2) / np.sum(np.array(inc_temp) < np.pi/2))
        print(slope_formated, len(inc_temp))  
        broj.append(len(inc_temp))
            
    except:
        pass
            
            

            
        