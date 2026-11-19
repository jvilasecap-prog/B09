import numpy as np
import matplotlib as plt

import geometry_properties as gp
import general_geometry as gg

#GET VALUES FROM AIRFOIL
filename = 'ms313.dat'
with open(filename) as f:
    title = f.readline().strip()
    n_upper, n_lower = (int(float(v)) for v in f.readline().split())

# loadtxt skips the blank lines automatically
data = np.loadtxt(filename, skiprows=2)   # shape (90, 2)

upper = data[:n_upper]                    # x, y of the upper surface
lower = data[n_upper:n_upper + n_lower]   # x, y of the lower surface

x_upper, y_upper = upper[:, 0], upper[:, 1]
x_lower, y_lower = lower[:, 0], lower[:, 1]

c_LE = 0.2
c_TE = 0.8 

yLE= np.array([np.interp(c_LE, x_lower, y_lower), np.interp(c_LE, x_upper, y_upper)])
yTE= np.array([np.interp(c_TE, x_lower, y_lower), np.interp(c_TE, x_upper, y_upper)])

#COMPUTE IXX DIAGRAMS

S = 26.65185874983523
AR = 11
b = np.sqrt(AR*S)
taper = 0.4
c_r = 2*S/(b*(1+taper))
c_t = taper*c_r

A = np.array([c_LE, yLE[1]])
B = np.array([c_TE, yTE[1]])
C = np.array([c_TE, yTE[0]])
D = np.array([c_LE, yLE[0]])\
t = 0.001#Preliminary value that can be change
A_array = np.array([])

Lines = gg.get_L_values(A, B, C, D, t)
Points = gg.get_A_values(A_array)
pos_cen = gg.global_Centroid(Lines, Points)
Ixx, Iyy = gg.global_I(Lines, Points, pos_cen)

#TRY TO PLOT THE IXX AND IYY
c = c_t + x*(c_r-c_t)/b
def c(x):
  c = c_t + x*(c_r-c_t)/b
plt.figure(figsize=(7, 4))
plt.plot(x, c_t + x*(c_r-c_t)/b, label=r"$f(x) = \sin(x)e^{-0.1x}$", linewidth=2)
plt.title("Damped Sine Wave")
plt.xlabel("x")
plt.ylabel("y")
plt.grid(True, linestyle="--", alpha=0.5)
plt.legend()
plt.tight_layout()
plt.show()
