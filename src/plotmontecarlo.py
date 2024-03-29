import numpy as np
import matplotlib.pyplot as plt

# Setting fonts for pretty plot
fonts = {
    "font.family": "serif",
    "axes.labelsize": 8,
    "font.size": 8,
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
}
plt.rcParams.update(fonts)

# Analytical integral value
analytical = 5 * np.pi ** 2 / 16 ** 2

# Loading data and defining errors
data_lambda_2 = np.loadtxt("montecarlo_lambda_2.txt", skiprows=1)
data_lambda_3 = np.loadtxt("montecarlo_lambda_3.txt", skiprows=1)
data_improved = np.loadtxt("montecarlo_improved.txt", skiprows=1)
dataO1 = np.loadtxt("montecarlo_paro1.txt", skiprows=1)
dataO2 = np.loadtxt("montecarlo_paro2.txt", skiprows=1)
dataO3 = np.loadtxt("montecarlo_paro3.txt", skiprows=1)


N_brute_lambda_2 = data_lambda_2[:, 0]
N_brute_lambda_3 = data_lambda_3[:, 0]
N_improved = data_improved[:, 0]
N_par = dataO1[:, 0]


error_brute_lambda_2 = np.abs((data_lambda_2[:, 1] - analytical) / analytical)
error_brute_lambda_3 = np.abs((data_lambda_3[:, 1] - analytical) / analytical)
error_improved = np.abs((data_improved[:, 1] - analytical) / analytical)

var_brute_lambda_2 = data_lambda_2[:, 2]
var_brute_lambda_3 = data_lambda_3[:, 2]
var_improved = data_improved[:, 2]

t_brute_lambda_2 = data_lambda_2[:, 3]
t_brute_lambda_3 = data_lambda_3[:, 3]
t_brute_par = data_lambda_3[:, 4]
t_improved = data_improved[:, 3]
t_improved_par = data_improved[:, 4]
t_O1 = dataO1[:, 3]
t_O2 = dataO2[:, 3]
t_O3 = dataO3[:, 3]

# Plotting errors
plt.tight_layout()
fig, ax = plt.subplots()
fig.set_size_inches(7.1014 / 2, 7.1014 / 2 * 3 / 4)
ax.loglog(N_brute_lambda_2, error_brute_lambda_2, label=r"Brute force, $\lambda = 2$")
ax.loglog(N_brute_lambda_3, error_brute_lambda_3, label=r"Brute force, $\lambda = 3$")
ax.loglog(N_improved, error_improved, label="Improved")
ax.grid()
ax.set_xlabel(r"$N$")
ax.set_ylabel("Relative Error")
plt.tight_layout()
ax.legend()
fig.savefig("../doc/Figures/error_monte_carlo.pdf", dpi=1000)

# Plotting variance
fig, ax = plt.subplots()
fig.set_size_inches(7.1014 / 2, 7.1014 / 2 * 3 / 4)
ax.loglog(N_brute_lambda_2, var_brute_lambda_2, label=r"Brute force, $\lambda = 2$")
ax.loglog(N_brute_lambda_3, var_brute_lambda_3, label=r"Brute force, $\lambda = 3$")
ax.loglog(N_improved, var_improved, label="Improved")
ax.grid()
ax.set_xlabel(r"$N$")
ax.set_ylabel("Variance")
plt.tight_layout()
ax.legend()

fig.savefig("../doc/Figures/variance_monte_carlo.pdf", dpi=1000)

# Plotting run time
fig, ax = plt.subplots()
fig.set_size_inches(7.1014 / 2, 7.1014 / 2 * 3 / 4)
ax.plot(N_brute_lambda_3, t_brute_lambda_3, label="Brute, 1 thread")
ax.plot(N_brute_lambda_3, t_brute_par, label="Brute, 2 threads")
ax.plot(N_improved, t_improved, label="Improved, 1 thread")
ax.plot(N_improved, t_improved_par, label="Improved, 2 threads")
ax.grid()
ax.set_xlabel(r"$N$")
ax.set_ylabel("CPU time [ms]")
plt.tight_layout()
ax.legend()
fig.savefig("../doc/Figures/cpu_time_monte_carlo.pdf", dpi=1000)

# Plotting run time for three compiler flags
fig, ax = plt.subplots()
fig.set_size_inches(7.1014 / 2, 7.1014 / 2 * 3 / 4)
ax.plot(N_par, t_O1, label="O1")
ax.plot(N_par, t_O2, label="O2")
ax.plot(N_par, t_O3, label="O3")
ax.grid()
ax.set_xlabel(r"$N$")
ax.set_ylabel("CPU time [ms]")
plt.tight_layout()
ax.legend()
fig.savefig("../doc/Figures/cpu_time_compilerflag.pdf", dpi=1000)
