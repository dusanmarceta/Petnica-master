import pandas as pd
import numpy as np
slope_values = np.arange(1, 4.2, 0.2)

time_shift = 184.00014287723752

for slope in slope_values:
    
    slope_formated = f"{round(slope, 1):.1f}"
    print(slope_formated)
    for i in range(1, 21):
        orbit_file = f'asteroids/A_M_{slope_formated}_1.0_orbit_file_{i}_.txt'
        
        # Učitaj fajl
        df = pd.read_csv(orbit_file)
        
        # Dodaj vrednost na poslednje dve kolone
        df.iloc[:, -2] += time_shift
        df.iloc[:, -1] += time_shift
        
        # Sačuvaj nazad u istom formatu (bez indexa, sa headerom, zarez kao separator)
        df.to_csv(orbit_file, index=False)