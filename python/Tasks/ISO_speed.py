import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.constants import GM_sun, au


v0 = 100e3
r0 = au.value

v0 = np.sqrt(2 * GM_sun.value / r0)

C = 0

r = np.linspace(1, 1000, 1000) * au.value

v = np.sqrt(2 * (C + GM_sun.value / r))

plt.plot(r/au.value, v/1000, 'k', linewidth = 4)

plt.xlabel('heliocentric distance (au)',fontsize = 48)
plt.ylabel('heliocentric speed (km/s)', fontsize = 48)
plt.xticks(fontsize = 40)
plt.yticks(fontsize = 40)
plt.grid()
plt.subplots_adjust(left=0.1, right=0.99, top=0.9, bottom=0.15)



            

            
        