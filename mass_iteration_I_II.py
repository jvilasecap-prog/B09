import pandas as pd
import math

#import python files as modules
import torenbeek_weight as wc2
import drag_computation2 as dc2
import ISA
#import the class II mass estimation

df = pd.read_excel('final_wing_planform.xlsx')
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
b_wing = wing['Wing span']
[P_cruise, T_cruise, rho_cruise] = ISA.calculate(41000, 0, 15, "ft")

tolerance = 5e-3
error = 1.0
#Iterate in a while loop until the MTOM error is less than 0.5% (arbitrary)
n=0
while (error > tolerance) and (n<100):
  n += 1
#Get the values for m_f, mtom, moem, C_D from an excel or the class II weight estimation
  CD_0 = dc2.drag_compute(mac_wing, wing["mac horizontal tail"], wing["mac vertical tail"], S_old, wing["area horizontal tail"], wing['area vertical tail'], wing["diameter engine"], 3, wing["Fuselage diameter"], wing["L1"], wing["L2"], wing["L3"], wing["Mach number cruise"], rho_cruise, T_cruise, wing["dynamic viscosity cruise"])
  oem = wc2.class_2_zero_fuel_weight(mtom_old, mzfw, wing["area horizontal tail"], b_wing, wing["Root chord"]*wing['t/c MS0313'], math.atan(-wing['Root chord']*(1-wing['Taper ratio'])/(2*b_wing)), 2, wing["area horizontal tail"], wing["area vertical tail"], wing['sweep H'], wing['sweep V'], 236.03, wing['H tail arm'], wing['Fuselage diameter'], wing['Fuselage diameter'], math.pi*(wing['Fuselage diameter']**2)/4)
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
  #compute the error from the optimal solution
  error = abs(mtom_new - mtom_old)/mtom_old 

  S_new = mtom_new*wing['Wing loading (initial, before any fuel has been burned)']
  b_wing = b_wing*math.sqrt(S_new/S_old)
  mac_wing = mac_wing*math.sqrt(S_new/S_old)

#the update of the excel was deleted since more testing is probably needed before the final reliable values are achieved
  mtom_old = mtom_new
  S_old = S_new

print(f"The final values are m_f = {m_f} kg\nm_tom_new = {mtom_new} kg\nm_oem = {oem} kg\nC_D0 = {CD_0}")
print(f"The final wing values are S = {S_new} m^2\n")


