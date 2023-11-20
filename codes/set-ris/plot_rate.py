import numpy as np

import matplotlib.pyplot as plt

from scipy.io import loadmat

# Load data
data = loadmat('set-ris.mat')

# Get probability vector
proba_vec = data['proba_vec'].squeeze()

# Get average delay
avg_delay = data['avg_delay']
avg_delay = avg_delay.mean(axis=1)

# Get average rate
rate = data['rate']
rate = rate.mean(axis=1)

# Lyapunov Constant
V1 = np.array([1e8, 1e9, 2e9, 3e10])
labels = ['1e8', '1e9', '2e9', '3e10']

fig, axes = plt.subplots(ncols=2)

for vv, V in enumerate(V1):
    axes[0].plot(proba_vec, avg_delay[:, vv], label='V =' + str(labels[vv]))
    axes[1].plot(proba_vec, rate[:, vv])

axes[0].set_xlabel('Probability of losing SET-R')
axes[0].set_ylabel('Average Delay')

axes[1].set_xlabel('Probability of losing SET-R')
axes[1].set_ylabel('Average Rate')

axes[1].set_yscale('log')

axes[0].legend()

plt.tight_layout()

plt.savefig('figs/set-ris_D500.pdf')

plt.show()

# # Number of users
# K = 2
#
# # Get rates
# rates = []
# for file in files:
#
#
#
#     rates.append(rate)
#
# rates = np.array(rates)
#

#
# #ax.plot(prob_vec,)
#
# # Access the data in Python
# #data = mat_data['data']
# #print(data)