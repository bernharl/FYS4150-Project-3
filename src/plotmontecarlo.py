import numpy as np
import matplotlib.pyplot as plt

# Setting fonts for pretty plot
fonts = {
    "font.family": "serif",
    "axes.labelsize": 10,
    "font.size": 10,
    "legend.fontsize": 10,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
}
plt.rcParams.update(fonts)

# Analytical integral value
analytical = 5 * np.pi**2 / 16**2

# Loading data and defining errors
data_brute = np.loadtxt("montecarlo.txt", skiprows = 1)
data_improved = np.loadtxt("montecarlo_improved.txt", skiprows = 1)
dataO1 = np.loadtxt("montecarlo_paro1.txt", skiprows = 1)
dataO2 = np.loadtxt("montecarlo_paro2.txt", skiprows = 1)
dataO3 = np.loadtxt("montecarlo_paro3.txt", skiprows = 1)


N_brute = data_brute[:, 0]
N_improved = data_improved[:, 0]
N_par = dataO1[:, 0]


error_brute = np.abs((data_brute[:, 1] - analytical) / analytical)
error_improved = np.abs((data_improved[:, 1] - analytical) / analytical)

var_brute = data_brute[:, 2]
var_improved = data_improved[:, 2]

t_brute = data_brute[:, 3]
t_improved = data_improved[:, 3]
t_O1 = dataO1[:, 3]
t_O2 = dataO2[:, 3]
t_O3 = dataO3[:, 3]

# Plotting errors
fig, ax = plt.subplots()
ax.loglog(N_brute, error_brute, label = "Brute force")
ax.loglog(N_improved, error_improved, label = "Improved")
ax.set_xlabel(r"$N$")
ax.set_ylabel("Relative Error")
ax.legend()
fig.savefig("../doc/Figures/error_monte_carlo.eps", dpi=1000)

# Plotting variance
fig, ax = plt.subplots()
ax.loglog(N_brute, var_brute, label = "Brute force")
ax.loglog(N_improved, var_improved, label = "Improved")
ax.set_xlabel(r"$N$")
ax.set_ylabel("Variance")
ax.legend()
fig.savefig("../doc/Figures/variance_monte_carlo.eps", dpi=1000)

# Plotting run time
fig, ax = plt.subplots()
ax.plot(N_brute, t_brute, label = "Brute force")
ax.plot(N_improved, t_improved, label = "Improved" )
ax.set_xlabel(r"$N$")
ax.set_ylabel("CPU time [ms]")
ax.legend()
fig.savefig("../doc/Figures/cpu_time_monte_carlo.eps", dpi=1000)

# Plotting run time for three compiler flags
fig, ax = plt.subplots()
ax.plot(N_par, t_O1, label = "O1")
ax.plot(N_par, t_O2, label = "O2" )
ax.plot(N_par, t_O3, label = "O3")
ax.set_xlabel(r"$N$")
ax.set_ylabel("CPU time [ms]")
ax.legend()
fig.savefig("../doc/Figures/cpu_time_compilerflag.eps", dpi=1000)
plt.show()