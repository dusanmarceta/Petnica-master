import numpy as np
import matplotlib.pyplot as plt
#from synthetic_population import synthetic_population_shell

sample = int(1e5)

# Sun
u_Sun=1e4    # u-component of the solar motion with respect to the LSR (m/s)
v_Sun=1.1e4  # w-component of the solar motion with respect to the LSR (m/s)
w_Sun=7e3    # w-component of the solar motion with respect to the LSR (m/s)

# ISOs
sigma_u=3.1e4 # dispersion of the u-component of ISOs with respect to the LSR when far from the Sun (m/s)
sigma_v=2.3e4 # dispersion of the v-component of ISOs with respect to the LSR when far from the Sun (m/s)
sigma_w=1.6e4 # dispersion of the w-component of ISOs with respect to the LSR when far from the Sun (m/s)


u = np.random.randn(sample) * sigma_u - u_Sun
v = np.random.randn(sample) * sigma_v - v_Sun
w = np.random.randn(sample) * sigma_w - w_Sun

V = np.sqrt(u**2 + v**2 + w**2)

plt.subplot(231)
plt.hist(u, 100, density = True)
plt.xlabel('u (m/s)', fontsize = 12)
plt.title('u - component', fontsize = 16, fontweight = 'bold', color = 'r')

plt.subplot(232)
plt.hist(v, 100, density = True)
plt.xlabel('v (m/s)', fontsize = 12)
plt.title('v - component', fontsize = 16, fontweight = 'bold', color = 'r')

plt.subplot(233)
plt.hist(w, 100, density = True)
plt.xlabel('w (m/s)', fontsize = 12)
plt.title('w - component', fontsize = 16, fontweight = 'bold', color = 'r')


plt.subplot(235)
plt.hist(V, 100, density = True)
plt.xlabel('speed (m/s)', fontsize = 12)
plt.title('speed', fontsize = 16, fontweight = 'bold', color = 'r')

plt.subplots_adjust(hspace=0.6)

print(f'90th percentile of speed is {np.round(np.percentile(V, 90), 2)} m/s')
print(f'99th percentile of speed is {np.round(np.percentile(V, 99), 2)} m/s')