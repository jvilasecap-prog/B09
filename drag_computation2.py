import math

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
def drag_compute(mac_wing, mac_H, mac_V, S_wing, S_H, S_V, D_engine, l_engine, D_fuselage, L1, L2, L3, M, rho_cruise, T_cruise, visc_cruise):
  #cf = 0.0045 #Light aircraft twin engine

#Computation of the wet surface
#The Surface area of the nacelle is negligible since the engine is attached to the fuselage
  S_wet_wing = 1.07*2*S_wing
  S_wet_V = 1.05*2*S_V
  S_wet_H = 1.05*2*S_H
  S_wet_engine = 2*math.pi*D_engine*l_engine
#Computation of wet surface of the fuselage
  term_bracket = (1 / (3 * L1**2)) * (((4 * L1**2 + (D_fuselage**2) / 4) ** 1.5) - ((D_fuselage**3) / 8))
  remaining_terms = -D_fuselage + 4 * L2 + 2 * math.sqrt(L3**2 + (D_fuselage**2) / 4)
    
  S_wet_fuselage = (math.pi * D_fuselage / 4) * (term_bracket + remaining_terms)
  #Computation of speed of sound
  a = math.sqrt(1.4*287*T_cruise)
  #k-value for smooth paint used for both wings and fuselage
  k = 0.634e-5
  #aircraft is assumed to be a classic production metal as a conserative estimation
  lam_part_wing = 0.1
  lam_part_tails = 0.1
  lam_part_fuselage = 0
  lam_part_engine = 0 #assumed 0 as a conservatime estimate
#computation of Re
  Re_wing = compute_Re(rho_cruise, M, a, mac_wing, visc_cruise, k)
  Re_H = compute_Re(rho_cruise, M, a, mac_H, visc_cruise, k)
  Re_V = compute_Re(rho_cruise, M, a, mac_V, visc_cruise, k)
  Re_fuselage = compute_Re(rho_cruise, M, a, L1+L2+L3, visc_cruise, k)
  Re_engine = compute_Re(rho_cruise, M, a, l_engine, visc_cruise, k)

#Computation of C_F
  C_f_wing = compute_Cf(lam_part_wing, Re_wing, M)
  C_f_H = compute_Cf(lam_part_tails, Re_H, M)
  C_f_V = compute_Cf(lam_part_tails, Re_V, M)
  C_f_fuselage = compute_Cf(lam_part_fuselage, Re_fuselage, M)
  C_f_engine = compute_Cf(lam_part_engine, Re_engine, M)

#computation of FF

  #For MS0313 the data 
  x_c_m_wing = 0.375 #location of maximum thickness
  t_c_wing = 0.131 #maximum thickness
  Lambda_m_wing = 0
  FF_wing = (1 + (0.6 / x_c_m_wing) * (t_c_wing) + 100 * (t_c_wing ** 4)) * (1.34 * (M ** 0.18) * (math.cos(Lambda_m_wing) ** 0.28))
  #For both tails a NACA0012 is assumed
  x_c_m_tail = 0.3
  t_c_tail = 0.12
  Lambda_m_H = 30*math.pi/180
  Lambda_m_V = 40*math.pi/180
  FF_H = (1 + (0.6 / x_c_m_tail) * (t_c_tail) + 100 * (t_c_tail ** 4)) * (1.34 * (M ** 0.18) * (math.cos(Lambda_m_H) ** 0.28))  
  FF_V = (1 + (0.6 / x_c_m_tail) * (t_c_tail) + 100 * (t_c_tail ** 4)) * (1.34 * (M ** 0.18) * (math.cos(Lambda_m_V) ** 0.28))  
  
  f_fuselage = (L1 + L2 + L3)/D_fuselage

  FF_fuselage = 1 + (60 / (f_fuselage ** 3)) + (f_fuselage / 400)
  f_engine = l_engine/D_engine
  FF_engine = 1 + 0.35/f_engine
  
  #IF computation
  #Due to T-tail
  IF_H = 1.04
  IF_V = 1.04

  IF_wing = 1
  #IF_fuselage is due to the pressence of the engines
  IF_fuselage = 1.3
  IF_engine = 1.3
  #wave drag considered negligible since M < M_DD

  C_D = (C_f_wing*FF_wing*IF_wing*S_wet_wing + C_f_H*FF_H*IF_H*S_wet_H + C_f_V*FF_V*IF_V*S_wet_V + C_f_fuselage*FF_fuselage*IF_fuselage*S_wet_fuselage + C_f_engine*FF_engine*IF_engine*S_wet_engine)/S_wing
  return C_D

  #to check the file
  #print(drag_compute(mac_wing, mac_H, mac_V, S_wing, S_H, S_V, D_engine, l_engine, D_fuselage, L1, L2, L3, M, rho_cruise, T_cruise, visc_cruise))