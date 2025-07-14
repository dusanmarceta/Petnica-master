import numpy as np
import matplotlib.pyplot as plt



'''
d_0:  reference ISO diameter for which n0 is defined (m)
d: array of diemetars for where power law for size frequency distribution (SFD) changes slope. This array also includes
   minimum and maximum diameter od the population (m). If this array is empty (default) the function does not calculate sizes of the objects 
   and takes n0 as the total number-density 
alpha: array of slopes of the SFD
    

Output (synthetic samples of orbital elements):
D_s: diameters of objects
'''


sample = int(1e6)
d_min = 10
d_max = 1e5
d_0 = 100
n_0 = 0.1
q = -2.5


n_max = n_0 * (d_min/d_0)**q
n_min = n_0 * (d_max/d_0)**q

n_total = n_max - n_min
x = np.random.random(sample) * n_total

D=d_min*((n_max-x)/n_max)**(1/q)


print(f'90th percentile of diameter is {np.round(np.percentile(D, 90), 2)} m')
print(f'99th percentile of diameter is {np.round(np.percentile(D, 99), 2)} m')