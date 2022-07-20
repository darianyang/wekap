"""
Plot the rates.
"""

import numpy as np
import matplotlib.pyplot as plt
import os

scheme = "../20-100conWE_lt36C2"
iterations = 500

os.chdir(scheme)

# load/gen data
raw_rates = np.load("rates.npy")
rates = np.load("raw_rates.npy")
xs = np.arange(1, iterations + 1, 1)

fig, ax = plt.subplots()
#ax.set_yscale('log')

ax.scatter(xs, raw_rates)

plt.show()
