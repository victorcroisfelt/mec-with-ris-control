import numpy as np

import matplotlib.pyplot as plt

from scipy.io import loadmat

# Load data
data = loadmat('data/set-ris.mat')

# Get probability vector
proba_vec = data['proba_vec']

# Get average delay
avg_delay = data['avg_delay']


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
# fig, ax = plt.subplots()
#
# ax.plot(prob_vec, rates[:, 0], label='UE 1')
# ax.plot(prob_vec, rates[:, 1], label='UE 2')
#
# ax.set_xlabel('probability of losing INI-R')
# ax.set_ylabel('UL rate')
#
# ax.legend()
#
# plt.tight_layout()
# plt.show()
#
# #ax.plot(prob_vec,)
#
# # Access the data in Python
# #data = mat_data['data']
# #print(data)