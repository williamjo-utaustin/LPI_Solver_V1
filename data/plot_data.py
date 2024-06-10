import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt("scratch.csv", delimiter = " ")
plt.semilogx(data[:,0], data[:,1])
plt.xlim(1E2, 1E8)
plt.ylim(0,0.10)

plt.show()

plt.loglog(data[:,0], data[:,2])
plt.loglog(data[:,0], data[:,3])
plt.xlim(1E2, 1E8)
plt.ylim(1E0, 1E-7)
plt.show()