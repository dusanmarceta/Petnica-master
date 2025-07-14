import numpy as np
import os

slope_values = np.arange(1, 4.2, 0.2)

for slope in slope_values:
    slope_formatted = f"{round(slope, 1):.1f}"
    for i in range(1, 21):
        pp_file = f'asteroids/A_M_{slope_formatted}_1.0_pp_file_{i}_.txt'

        if os.path.exists(pp_file):
            # Učitaj sve linije
            with open(pp_file, 'r') as f:
                lines = f.readlines()

            # Zameni prvi red
            lines[0] = 'ObjID,H_r,u-r,g-r,i-r,z-r,y-r,GS\n'

            # Upisi nazad
            with open(pp_file, 'w') as f:
                f.writelines(lines)

            print(f"Updated: {pp_file}")
        else:
            print(f"File not found: {pp_file}")
