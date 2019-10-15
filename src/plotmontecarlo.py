import numpy as np 
import matplotlib.pyplot as plt 

analytical = 5 * np.pi**2 / 16**2

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


fig, ax = plt.subplots()
ax.loglog(N_brute, error_brute, label = "Brute")
ax.loglog(N_improved, error_improved, label = "Improved")
ax.set_xlabel(r"$N$")
ax.set_ylabel("relative error")
ax.legend()

fig, ax = plt.subplots()
ax.loglog(N_brute, var_brute, label = "Brute")
ax.loglog(N_improved, var_improved, label = "Improved")
ax.set_xlabel(r"$N$")
ax.set_ylabel("Variance")
ax.legend()

fig, ax = plt.subplots()
ax.plot(N_brute, t_brute, label = "Brute")
ax.plot(N_improved, t_improved, label = "Improved" )
ax.set_xlabel(r"$N$")
ax.set_ylabel("CPU time [ms]")
ax.legend()


fig, ax = plt.subplots()
ax.plot(N_par, t_O1, label = "O1")
ax.plot(N_par, t_O2, label = "O2" )
ax.plot(N_par, t_O3, label = "O3")
ax.set_xlabel(r"$N$")
ax.set_ylabel("CPU time [ms]")
ax.legend()
plt.show()