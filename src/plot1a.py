import numpy as np
import matplotlib.pyplot as plt

# Setting fonts for pretty plot
fonts = {
    "font.family": "serif",
    "axes.labelsize": 16,
    "font.size": 16,
    "legend.fontsize": 16,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
}
plt.rcParams.update(fonts)

# Reading data, defining errors

data = np.loadtxt("Exercise_a_b.txt", skiprows = 1)
N = data[:, 0]
error_gauleg = data[:, 1]
error_gauss_laguerre = data[:, 2]

fig, ax = plt.subplots()
# The golden ratio for nice plots size, also based on width of tex document
# fig.set_size_inches(2 * 2.9, 2 * 1.81134774961)

# Plotting data
ax.semilogy(N, error_gauleg, label = "Brute force")
ax.semilogy(N, error_gauss_laguerre, label = "Improved")
ax.semilogy(N, np.ones_like(N) * 1e-2, label = "Three digit precision")
ax.grid()
ax.set_xlabel(r"$N$")
ax.set_ylabel(r"Absolute error")
fig.tight_layout()
ax.legend(loc=1)
# Saving high quality figure
fig.savefig("../doc/Figures/exercise_a_b.pdf", dpi=1000)