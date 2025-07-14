import numpy as np
import matplotlib.pyplot as plt
from astropy.constants import au

V_cut = 24.5
D_max = 60000
albedo_max = 0.1
v_max = 100000
survey_time = 10
n_max = 300

H = 15.618 - 2.5 * np.log10(albedo_max) - 5 * np.log10(D_max / 1000)

C = 10 ** ((V_cut - H) / 5)
distance = (1 + np.sqrt(1 + 4 * C)) / 2

print(np.shape(distance))



r_model = distance + (survey_time * 365.25*86400 * v_max)/au.value


N = 4/3 * r_model**3 * np.pi * n_max


# how r_model depands on D_max and v_max



D = np.linspace(10, 60000, 100)
v = np.linspace(100, 100000, 50)

D, v = np.meshgrid(D, v)


H = 15.618 - 2.5 * np.log10(albedo_max) - 5 * np.log10(D / 1000)

C = 10 ** ((V_cut - H) / 5)
distance = (1 + np.sqrt(1 + 4 * C)) / 2

r_model_array = distance + (survey_time * 365.25*86400 * v)/au.value


plt.figure()


contour = plt.contourf(D/1000, v/1000, r_model_array, 100, cmap='jet')
plt.xlabel('D (km)', fontsize=40)
plt.ylabel('v (km/s)', fontsize=40)

plt.xticks(fontsize=36)
plt.yticks(fontsize=36)

cbar = plt.colorbar(contour)
cbar.set_label('model sphere radius (au)', fontsize=40)
cbar.ax.tick_params(labelsize=36)  # OVA LINIJA povećava cbar tickove
plt.grid()

plt.subplots_adjust(left=0.11, right=0.99, top=0.98, bottom=0.12)




plt.figure()
contour = plt.contourf(D/1000, v/1000, r_model_array, 100, cmap='jet')
plt.xlabel('D (km)', fontsize=40)
plt.ylabel('v (km/s)', fontsize=40)

plt.xticks(fontsize=36)
plt.yticks(fontsize=36)

cbar = plt.colorbar(contour)
cbar.set_label('model sphere radius (au)', fontsize=40)
cbar.ax.tick_params(labelsize=36)  # OVA LINIJA povećava cbar tickove
plt.grid()

plt.subplots_adjust(left=0.11, right=0.99, top=0.98, bottom=0.12)



