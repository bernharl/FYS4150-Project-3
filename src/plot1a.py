import numpy as np 
import matplotlib.pyplot as plt 


data = np.loadtxt("Exercise_a_b.txt", skiprows = 1)
N = data[:, 0]
error = data[:, 1]
plt.plot(N, error)
plt.show()