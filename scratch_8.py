import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append('/Users/williamjo/Documents/LPI/Codes/plume_regolith_solver_v1/src')


import var_constants as cs

u_ej = np.linspace(0, 8000, 100)


v_ej_bar = (u_ej**2) / (cs.r_moon * cs.g)
r_max_1 = 2 * (cs.r_moon) * np.arctan((v_ej_bar**2 * np.sin(3*np.pi/180) * np.cos(3 * np.pi/180))/(1 - v_ej_bar**2 * (np.cos(3*np.pi/180))**2))


r_max_2 = (u_ej**2 * np.sin(2 * 3 * np.pi/180))/ (cs.g)

print(5442068.267824446*2)

for i in range(0,100):
    print(u_ej[i], r_max_1[i], r_max_2[i], 10920176)


vesc = np.sqrt(2 * cs.g * cs.r_moon)
print(vesc)
print(2 * np.pi * cs.r_moon)
plt.semilogy(u_ej, r_max_1)
plt.semilogy(u_ej, r_max_2)
plt.show()



