import numpy as np

import matplotlib.pyplot as plt

from scipy.io import loadmat

# Files
files = ['100_k1', '125_k1', '150_k1', '200_k1', '250_k1', '300_k1']
#files = ['baseline', 'ini-ris', 'set-ris', 'set-ue']

labels = ['Baseline', 'INI-RIS', 'SET-RIS', 'INI-UE', 'SET-UE']
#labels = ['Baseline', 'INI-RIS', 'SET-RIS', 'SET-UE']

# Constant
V1 = np.array([1e8, 1e9, 2e9, 3e10])

fig, ax = plt.subplots()

y_axis = []
x_axis = [100, 125, 150, 200, 250, 300]

# Go through the files
for ff, file in enumerate(files):

    # Load data
    data = loadmat(file + '.mat')

    # Compute average and take maximum
    avg_delay = np.mean(data['avg_delay'], axis=0).min()

    y_axis.append(avg_delay)

ax.plot(x_axis, y_axis)


#     if ff == 0:
#
#         # Get average delay
#         avg_delay = data['avg_delay']
#         avg_delay = avg_delay.mean(axis=0)
#
#         ax.plot(np.linspace(0, 1, 11), np.max(avg_delay/1000, axis=-1) * np.ones(11), color='black', label=labels[ff])
#
#     else:
#
#         # Get probability vector
#         proba_vec = data['proba_vec'].squeeze()
#
#         # Get average delay
#         avg_delay = data['avg_delay']
#         avg_delay = avg_delay.mean(axis=1)
#
#         if labels[ff] == 'INI-UE':
#             ax.plot(proba_vec, np.flip(np.max(avg_delay/1000, axis=-1)), label=labels[ff])
#         else:
#             ax.plot(proba_vec, (np.max(avg_delay / 1000, axis=-1)), label=labels[ff])
#
#
ax.set_yscale('log')

ax.set_xlabel(r'Slot duration, $\tau$')
ax.set_ylabel('Average delay [ms]')

ax.legend()

plt.tight_layout()
plt.show()
