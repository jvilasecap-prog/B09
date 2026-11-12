import pandas as pd
import math

#import python files as modules
import torenbeek_weight as wc2
import drag_computation2 as dc2
import ISA
import drag_polar
#import the class II mass estimation

df = pd.read_excel('final_wing_planform_corrected.xlsx')
wing_unclean = df.set_index('Parameter')['value'].to_dict()
wing = {}
for key, value in wing_unclean.items():
    try:
        wing[key] = float(value)
    except (ValueError, TypeError): # for strings or blanks
        wing[key] = value 
#print(wing)

# mass estimate from Class I:


mtom_old = wing['MTOM']
oem = wing['OEM']
m_f = wing['Fuel mass']
mzfw = oem + wing['Maximum payload']

S_old = wing['Wing area']
mac_wing = wing['mac wing']
c_r_wing = wing['Root chord']
b_wing = wing['Wing span']

S_H = wing['area horizontal tail']
mac_H = wing['mac horizontal tail']
S_V = wing['area vertical tail']
mac_V = wing['mac vertical tail']

[P_cruise, T_cruise, rho_cruise] = ISA.calculate(41000, 0, 15, "ft")

tolerance = 5e-4
error = 1.0
#Iterate in a while loop until the MTOM error is less than 0.05% (arbitrary)
n=0
wing_volume = drag_polar.calculate_volume("MS(1)-0313.txt", wing["Aspect ratio"], S_old, wing["Taper ratio"], 0.2, 0.8)
fuel_volume = m_f/805 
print(f"The final values are m_f = {m_f} kg\nmtom_new = {mtom_old} kg\noem = {oem} kg")
print(f"The final wing values are S = {S_old} m^2\n")
print(f"The available volume of the wings is {wing_volume} and the requried is {fuel_volume}")
print(f"The drag coefficient{dc2.drag_compute(
    mac_wing, 
    mac_H, 
    mac_V, 
    S_old, 
    S_H, 
    S_V, 
    wing["diameter engine"], 
    3, 
    wing["Fuselage diameter"], 
    wing["L1"], 
    wing["L2"], 
    wing["L3"], 
    wing["Mach number cruise"], 
    rho_cruise, 
    T_cruise, 
    wing["dynamic viscosity cruise"]
  )}")

print(wc2.w_wing(
    mzfw,
    S_old,
    b_wing,
    c_r_wing*wing['t/c MS0313'],
    math.atan(-c_r_wing*(1-wing['Taper ratio'])/(2*b_wing)),
    4,
))
print(wc2.w_fuselage(
    230,
    wing['H tail arm'],
    wing['Fuselage diameter'],
    wing['Fuselage diameter'],
     math.pi*wing['Fuselage diameter']*wing['Fuselage length'],
))
print(wc2.w_control_surfaces(mtom_old))
print(wc2.w_landing_gear(mtom_old))

a=0.947
#dummy variables to start the while loop
#wing_volume = 0
#fuel_volume = 1
#while (wing_volume < fuel_volume)and (n<100):
  #a -= 0.005
W_S = a*wing['Wing loading (initial, before any fuel has been burned)']
while (error > tolerance) and (n<100):
  n += 1
  #Get the values for m_f, mtom, moem, C_D from an excel or the class II weight estimation
  CD_0 = dc2.drag_compute(
    mac_wing, 
    mac_H, 
    mac_V, 
    S_old, 
    S_H, 
    S_V, 
    wing["diameter engine"], 
    1.923, 
    wing["Fuselage diameter"], 
    wing["L1"], 
    wing["L2"], 
    wing["L3"], 
    wing["Mach number cruise"], 
    rho_cruise, 
    T_cruise, 
    wing["dynamic viscosity cruise"]
  )
  oem = wc2.class_2_operating_empty_weight(mtom_old, mzfw, S_old, b_wing, c_r_wing*wing['t/c MS0313'], math.atan(-c_r_wing*(1-wing['Taper ratio'])/(2*b_wing)), wing['Ultimate load factor'], S_H, S_V, wing['sweep H'], wing['sweep V'], 236.03, wing['H tail arm'], wing['Fuselage diameter'], wing['Fuselage diameter'], math.pi*wing['Fuselage diameter']*wing['Fuselage length'], 5.08571194640947, 4.100936556)
  #CHECK EXCEL VALUES WITH GROUP
  mzfw = oem + wing['Maximum payload']
  #The m_f can be obtained with the range of the aircraft, however this may give issues in the convergence of the code since the aircraft is becoming heavier while the dimensional components are not growing
  #For that reason reitaration of the aircaft dimensional variables is probably needed
  #The fuel of the aircraft has to fit in the wings, there is a function to get the wing volume in a previous file

  #Computation of m_f
  e_p = 43*10**6 # [J/kg]
  bypass_ratio = 5.5 # (ADSEE)
  psi = wing["Psi"] # (ADSEE)
  phi = wing["Phi"] # (ADSEE)
  AR = wing["Aspect ratio"] # design choice

  e = 1/(math.pi*AR*psi+ 1/phi)
  TSFC = 22*bypass_ratio**(-0.19)

  SoS = math.sqrt(1.4 * ISA.R * T_cruise)
  V_CR = wing["Mach number cruise"] * SoS

  jet_efficiency = (V_CR/(TSFC*10**-6))/e_p
  L_D_max_cruise = 0.5* math.sqrt((math.pi*AR*e)/(CD_0))

  R_eq_res = 805165.68 # [m] (constant adsee calc based on fixed TLARs)
  R_lost = 1/0.7*L_D_max_cruise*(12496.8+(V_CR**2/(2*9.80665)))
  R_eq = (6100*1000 + R_lost)*1.05 + R_eq_res

  m_f_mtom_ratio = 1 - math.exp((-R_eq*9.80665)/(jet_efficiency*L_D_max_cruise*e_p))
  #computation of m_mtom_new
  mtom_new = mzfw/(1-m_f_mtom_ratio)

  m_f = mtom_new*m_f_mtom_ratio
  m_oem_mtow_ratio = oem/mtom_new

  #compute the error from the optimal solution
  error = abs(mtom_new - mtom_old)/mtom_old 
  #update of the wing values
  S_new = mtom_new*9.80665/W_S
  b_wing = b_wing*math.sqrt(S_new/S_old)
  mac_wing = mac_wing*math.sqrt(S_new/S_old)
  c_r_wing = c_r_wing*math.sqrt(S_new/S_old)
  S_V = S_V*(S_new/S_old)**(3/2) #The tail volume scales both with the wing area and the wingspan or chord
  S_H = S_H*(S_new/S_old)**(3/2)
  mac_V = mac_V*math.sqrt(S_new/S_old)
  mac_H = mac_H*math.sqrt(S_new/S_old)

  mtom_old = mtom_new
  S_old = S_new

  wing_volume = drag_polar.calculate_volume("MS(1)-0313.txt", wing["Aspect ratio"], S_new, wing["Taper ratio"], 0.2, 0.8)
  fuel_volume = m_f/805




#the update of the excel was deleted since more testing is probably needed before the final reliable values are achieved

#CHECK THAT THE FUEL CAN BE STORED IN THE WINGS
print(f"The final values are m_f = {m_f} kg\nmtom_new = {mtom_new} kg\noem = {oem} kg\nC_D0 = {CD_0}")
print(f"The final wing values are S = {S_new} m^2\n")
print(f"The available volume of the wings is {wing_volume} and the requried is {fuel_volume}")
print(f"Ratio OEM/MTOM{oem/mtom_new}")
print(f"W/S{W_S, a}")

