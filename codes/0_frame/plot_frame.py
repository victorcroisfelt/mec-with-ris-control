import numpy as np

import matplotlib.pyplot as plt

from scipy.io import loadmat

from scipy.interpolate import UnivariateSpline

fig, ax = plt.subplots()

for rr in range(2):

    if rr == 0:
        # Files
        files = ['100_k1', '125_k1', '150_k1', '200_k1', '250_k1', '300_k1']
        #files = ['baseline', 'ini-ris', 'set-ris', 'set-ue']

        labels = ['Baseline', 'INI-RIS', 'SET-RIS', 'INI-UE', 'SET-UE']
        #labels = ['Baseline', 'INI-RIS', 'SET-RIS', 'SET-UE']

        # Constant
        V1 = np.array([1e8, 1e9, 2e9, 3e10])

        y_axis = []
        x_axis = [100, 125, 150, 200, 250, 300]

        # Go through the files
        for ff, file in enumerate(files):

            # Load data
            data = loadmat(file + '.mat')

            # Compute average and take maximum
            avg_delay = np.mean(data['avg_delay'], axis=0).min()

            y_axis.append(avg_delay)

        # Create a cubic spline interpolation
        f = UnivariateSpline(x_axis, y_axis, k=1)

        # Define a finer grid for plotting
        x_interp = np.linspace(100, 300, 100)
        y_interp = f(x_interp)

        ax.plot(x_interp, y_interp, label='K=1')

    else:
        # Files
        files = ['200_k4', '300_k4', '400_k4']
        #files = ['baseline', 'ini-ris', 'set-ris', 'set-ue']

        #labels = ['Baseline', 'INI-RIS', 'SET-RIS', 'INI-UE', 'SET-UE']
        #labels = ['Baseline', 'INI-RIS', 'SET-RIS', 'SET-UE']

        # Constant
        V1 = np.array([1e8, 1e9, 2e9, 3e10])



        y_axis = []
        x_axis = [200, 300, 400]

        # Go through the files
        for ff, file in enumerate(files):

            # Load data
            data = loadmat(file + '.mat')

            # Compute average and take maximum
            avg_delay = np.mean(data['avg_delay'], axis=0).min()

            y_axis.append(avg_delay)

        # Create a cubic spline interpolation
        f = UnivariateSpline(x_axis, y_axis, k=1)

        # Define a finer grid for plotting
        x_interp = np.linspace(200, 400, 100)
        y_interp = f(x_interp)

        ax.plot(x_interp, y_interp, label='K=4')

ax.set_yscale('log')

ax.set_xlabel(r'Slot duration, $\tau$ [ms]')
ax.set_ylabel('Average delay [ms]')

ax.legend()

plt.tight_layout()

# import tikzplotlib
#
# tikzplotlib.save("frame.tex")

plt.show()
