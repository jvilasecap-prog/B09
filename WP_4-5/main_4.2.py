import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

import geometry_properties as gp
import general_geometry as gg

def main():
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
  t = 0.002#Preliminary value that can be changed

  n_stations = 50
  z_span = np.linspace(0, b / 2, n_stations)
  
  #M, T and E G values to test
  T = np.linspace(0, 500, n_stations)
  M = np.linspace(0, -250, n_stations)
  G = 1e9
  E = 2*G

  Ixx_list = np.zeros(n_stations)
  Iyy_list = np.zeros(n_stations)

  J_list = np.zeros(n_stations)
  i = 0
  for y in z_span: 
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

    


  v = gg.get_deflection(M, Ixx_list, E, z_span)
  theta = gg.get_twist(T, J_list, G, z_span)

  #TRY TO PLOT THE IXX AND IYY
  fig, ax = plt.subplots(figsize=(8, 5))
  ax.plot(z_span, Ixx_list, label=r"$I_{xx}$", linewidth=2)
  ax.plot(z_span, Iyy_list, label=r"$I_{yy}$", linewidth=2)
  ax.plot(z_span, J_list, label=r"$J$", linewidth=2)
  ax.set_xlabel("Spanwise position y [m]")
  ax.set_ylabel(r"Second moment of area [m$^4$]")
  ax.set_title(r"$I_{xx}$, $I_{yy}$ and $J$ along the wingspan")
  ax.grid(True, alpha=0.3)
  ax.legend()
  plt.tight_layout()
  plt.show()

  fig, ax = plt.subplots(figsize=(8, 5))
  ax.plot(z_span, theta, label=r"$theta(rad)$", linewidth=2)
  ax.plot(z_span, v, label=r"$deflection(m)$", linewidth=2)
  ax.set_xlabel("Spanwise position y [m]")
  ax.set_ylabel(r"Deflections")
  ax.set_title(r"twist and bending distribution along the wingspan")
  ax.grid(True, alpha=0.3)
  ax.legend()
  plt.tight_layout()
  plt.show()
if __name__ == "__main__":
    main()