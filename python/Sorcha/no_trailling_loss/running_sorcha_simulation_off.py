import os
import numpy as np
import os.path
import time

import sys
sys.path.append('..')


#from synthetic_population import synthetic_population
from auxiliary_functions import run_command

slope_values = np.arange(2, 4.4, 0.4)

cadence_file = '../cadence_file/baseline_v4.3.5_10yrs.db'


for slope in slope_values:
    slope_formated = f"{round(slope, 1):.1f}"
    for i in range(1, 21):
        orbit_file = '../asteroids/A_M_' + slope_formated + '_1.0_orbit_file_' + str(i) + '_.txt'
        pp_file = '../asteroids/A_M_' + slope_formated + '_1.0_pp_file_' + str(i) + '_.txt'
        output_file = '../results/no_trailling/A_M_' + slope_formated + '_1.0_output_file_' + str(i)
        
#        print('aaaaa')
#        command = f'nohup sorcha run -c css_configuration.ini -p {pp_file} -ob {orbit_file} -pd {cadence_file} -o ./ -t {output_file} > /dev/null 2>&1 &'
        

        command = f"nohup sorcha run -f -c ../Rubin_circular_approximation_off.ini -p {pp_file} -ob {orbit_file} -pd {cadence_file} -o ./ -t {output_file} > log_file.log 2>&1 &"
        
#        command = f"nohup sorcha run -c css_configuration.ini -p {pp_file} -ob {orbit_file} -pd {cadence_file} -o ./ -t {output_file} &"

        
        run_command(command)
        
        # waiting for Sorcha to finish the simulation
        while not os.path.exists(output_file + '.csv'):
            time.sleep(10)
            
        np.savetxt('progress.log', [f'slope = {slope_formated}, file number = {i}'], fmt='%s')