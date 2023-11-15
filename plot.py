import numpy as np

import matplotlib.pyplot as plt

from scipy.io import loadmat

import tikzplotlib

# Files
files = ['baseline', 'ini-ris2', 'set-ris2', 'ini-ue4', 'set-ue3']
labels = ['Baseline', 'INI-RIS', 'SET-RIS', 'INI-UE', 'SET-UE']

# Constant
V1 = np.array([1e8, 1e9, 2e9, 3e10])

fig, ax = plt.subplots()

# Go through the files
for ff, file in enumerate(files):

    # Load data
    data = loadmat('data/' + file + '.mat')

    if ff == 0:

        # Get average delay
        avg_delay = data['avg_delay']
        avg_delay = avg_delay.mean(axis=0)

        ax.plot(np.linspace(0, 1, 11), np.min(avg_delay/1000, axis=-1) * np.ones(11), color='black', label=labels[ff])

    else:

        # Get probability vector
        proba_vec = data['proba_vec'].squeeze()

        # Get average delay
        avg_delay = data['avg_delay']
        avg_delay = avg_delay.mean(axis=1)

        if labels[ff] == 'INI-UE':
            ax.plot(proba_vec, np.min(avg_delay/1000, axis=-1), label=labels[ff])
        else:
            ax.plot(proba_vec, (np.min(avg_delay / 1000, axis=-1)), label=labels[ff])


ax.set_xscale('log')

ax.set_xlabel('Probability of losing packet')
ax.set_ylabel('Average Delay [ms]')

ax.legend()

plt.tight_layout()

tikzplotlib.save("error_prob.tex")

plt.show()
