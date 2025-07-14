import numpy as np
import matplotlib.pyplot as plt
'''
Takes albedo and/or diameter data from WISE and generates a sample of arbitrary
size that imitates the WISE data
'''

import pandas as pd
from auxiliary_functions import imitate_sample



file='neowise_data/neowise_mainbelt.csv'
n_bins=20 # number of bins for the available data
N = 100000 # number of output objects

H, D, aV, aIR = (pd.read_csv(file, usecols=[3,11,13,15])).T.to_numpy()


variable=aV

# -0.999 is a flag for not having albedo
index = np.argwhere(variable==-0.999)
variable = np.delete(variable, index)


synthetic_sample=imitate_sample(variable, n_bins, N)
plt.hist(synthetic_sample, 100, density = True)
plt.xlabel('albedo', fontsize = 24)
print(f'90th percentile of albedo is {np.round(np.percentile(synthetic_sample, 90), 2)}')
print(f'99th percentile of albedo is {np.round(np.percentile(synthetic_sample, 99), 2)}')