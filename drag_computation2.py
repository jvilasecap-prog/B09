import math
import numpy as np
import pandas as pd

def compute_Re(rho, M, a, l, visc, k):
  #The aircraft is considered as high subsonic aircraft
  Re1 = rho*M*a*l/visc
  Re2 = 38.21*(l/k)**(1.053)
  return min(Re1, Re2)
def compute_Cf(lam_part, Re, M):
  C_f_lam = 1.328/math.sqrt(Re)
  C_f_tur = 0.455 / ((math.log10(Re) ** 2.58) * ((1 + 0.144 * (M ** 2)) ** 0.65))
  C_f = C_f_lam*lam_part + C_f_tur * (1-lam_part)
  return C_f
def drag_compute(mac_wing, mac_H, mac_V, S_wing, S_H, S_V, D_fuselage, L1, L2, L3, M, rho_cruise, T_cruise, visc_cruise):
  #cf = 0.0045 #Light aircraft twin engine

#Cmomputation of wet surface
  S_wet_wing = 1.07*2*S_wing
  S_wet_V = 1.05*2*S_V
  S_wet_H = 1.05*2*S_H

#Computation of wet surface of the fuselage
  term_bracket = (1 / (3 * L1**2)) * (((4 * L1**2 + (D_fuselage**2) / 4) ** 1.5) - ((D**3) / 8))
  remaining_terms = -D_fuselage + 4 * L2 + 2 * math.sqrt(L3**2 + (D_fuselage**2) / 4)
    
  S_w_fus = (math.pi * D_fuselage / 4) * (term_bracket + remaining_terms)
  #Computation of speed of sound
  a = math.sqrt(1.4*287*T_cruise)
  #k-value for smooth paint used for both wings and fuselage
  k = 0.634e-5
  #The aircraft is catalogued as a general aviation smooth molded composite
  lam_fus = 0.25
  lam_wing = 0.5
  #aircraft is assumed to be a classic production metal as a conserative estimation
  lam_part_wing = 0.1
  lam_part_tails = 0.1
  lam_part_fuselage = 0
#computation of Re
  Re_wing = compute_Re(rho_cruise, M, a, mac_wing, visc_cruise, k)
  Re_H = compute_Re(rho_cruise, M, a, mac_H, visc_cruise, k)
  Re_V = compute_Re(rho_cruise, M, a, mac_V, visc_cruise, k)
  Re_fuselage = compute_Re(rho_cruise, M, a, L1+L2+L3, visc_cruise, k)

#Computation of C_F
  C_f_wing = compute_Cf(lam_part_wing, Re_wing, M)
  C_f_H = compute_Cf(lam_part_tails, Re_H, M)
  C_f_V = compute_Cf(lam_part_tails, Re_V, M)
  C_f_fuselage = compute_Cf(lam_part_fuselage, Re_fuselage, M)

#computation of FF

  #For MS0313 the data 
  x_c_m_wing = 0.375 #location of maximum thickness
  t_c_wing = 0.131 #maximum thickness
  Lambda_m_wing = 0
  FF_wing = (1 + (0.6 / x_c_m_wing) * (t_c_wing) + 100 * (t_c_wing ** 4)) * (1.34 * (M ** 0.18) * (math.cos(Lambda_m_wing) ** 0.28))
  #For both tails a NACA0012 is assumed
  x_c_m_tail = 0.3
  t_c_tail = 0.12
  FF_H = (1 + (0.6 / x_c_m_tail) * (t_c_tail) + 100 * (t_c_tail ** 4)) * (1.34 * (M ** 0.18) * (math.cos(Lambda_m) ** 0.28))  
  
  f = (L1 + L2 + L3)/D_fuselage
  FF_fuselage = 1 + (60 / (f ** 3)) + (f / 400)