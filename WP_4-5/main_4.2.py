import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

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

#Points for unit chord with coordinate system LE
A = np.array([c_LE, yLE[1]])
B = np.array([c_TE, yTE[1]])
C = np.array([c_TE, yTE[0]])
D = np.array([c_LE, yLE[0]])

# new points (origin at B, x left, y down)
A_new = B - A
B_new = B - B   # -> [0, 0]
C_new = B - C
D_new = B - D
print(f"A {A} B {B} C {C} D {D}")
print(f"A {A_new} B {B_new} C {C_new} D {D_new}")
t = 0.001#Preliminary value that can be changed

n_stations = 50
y_span = np.linspace(0, b / 2, n_stations)

Ixx_list = np.zeros(n_stations)
Iyy_list = np.zeros(n_stations)

J_list = np.zeros(n_stations)
i = 0
for y in y_span: 
  #Scale the points
  c_local = c_r + (c_t - c_r) * (y / (b / 2))
  A = np.array([c_LE, yLE[1]])*c_local
  B = np.array([c_TE, yTE[1]])*c_local
  C = np.array([c_TE, yTE[0]])*c_local
  D = np.array([c_LE, yLE[0]])*c_local

  #Fill in A_array with something
  A_array = np.array([])

  Lines = gg.get_L_values(A, B, C, D, t)
  Points = gg.get_A_values(A_array)

  pos_cen = gg.global_Centroid(Lines, Points)
  Ixx, Iyy = gg.global_I(Lines, Points, pos_cen)
  J = gg.global_J(Lines, A, B, C, D)

  Ixx_list[i] = Ixx
  Iyy_list[i] = Iyy
  J_list[i] = J

  i += 1


  



#TRY TO PLOT THE IXX AND IYY
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(y_span, Ixx_list, label=r"$I_{xx}$", linewidth=2)
ax.plot(y_span, Iyy_list, label=r"$I_{yy}$", linewidth=2)
ax.set_xlabel("Spanwise position y [m]")
ax.set_ylabel(r"Second moment of area [m$^4$]")
ax.set_title(r"$I_{xx}$ and $I_{yy}$ along the wingspan")
ax.grid(True, alpha=0.3)
ax.legend()
plt.tight_layout()
plt.show()