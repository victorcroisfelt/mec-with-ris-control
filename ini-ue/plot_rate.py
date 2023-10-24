import numpy as np

import matplotlib.pyplot as plt

from scipy.io import loadmat

# Load data
data = loadmat('data/set-ris_500_45deg.mat')

# Get probability vector
proba_vec = data['proba_vec'].squeeze()

# Get average delay
avg_delay = data['avg_delay']
avg_delay = avg_delay.mean(axis=0)

# Constant
V1 = np.array([1e8, 1e9, 2e9, 3e10])

fig, ax = plt.subplots()

for vv, V in enumerate(V1):
    ax.plot(proba_vec, avg_delay[:, vv], label='V =' +str(V1[vv]))

ax.set_xlabel('Probability of losing SET-R')
ax.set_ylabel('Average Delay')

ax.legend()

plt.tight_layout()
plt.show()



breakpoint()
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