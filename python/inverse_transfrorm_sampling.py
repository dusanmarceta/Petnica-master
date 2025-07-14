import numpy as np
import matplotlib.pyplot as plt

sample_size = int(1e6)

ys = np.random.random(sample_size) 

xs = np.sqrt(ys) * np.pi

plt.figure()

plt.hist(xs, 100, density = True)
plt.title('Linear distribution', fontsize = 20)
plt.grid()

plt.figure()

xs = np.arccos(1-2*ys)
plt.hist(xs, 100, density = True)
plt.title('Sinusoidal distribution', fontsize = 20)
plt.grid()







