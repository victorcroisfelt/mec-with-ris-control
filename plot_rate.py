import numpy as np

import matplotlib.pyplot as plt

from scipy.io import loadmat

# Files
files = ['ini-ris', 'set-ris']
labels = ['INI-RIS', 'SET-RIS']

# Constant
V1 = np.array([1e8, 1e9, 2e9, 3e10])

fig, ax = plt.subplots()

# Go through the files
for ff, file in enumerate(files):

    # Load data
    data = loadmat('mat/' + file + '.mat')

    # Get probability vector
    proba_vec = data['proba_vec'].squeeze()

    # Get average delay
    avg_delay = data['avg_delay']
    avg_delay = avg_delay.mean(axis=1)

    ax.plot(proba_vec, np.max(avg_delay/1000, axis=-1), label=labels[ff])

#ax.set_xscale('log')

ax.set_xlabel('Probability of losing packet')
ax.set_ylabel('Average Delay [ms]')

ax.legend()

plt.tight_layout()
plt.show()
