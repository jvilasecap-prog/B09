import pandas as pd
import math

#import python files as modules
import torenbeek_weight
import drag_computation2
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


m_tom_old = wing['MTOM']
m_oem_old = wing['OEM']
m_f_old = wing['Fuel mass']
CD_0_old = wing['CD_0']
mzfw_old = m_oem_old + wing['Maximum payload']
[P_cruise, T_cruise, rho_cruise] = ISA.calculate(41000, 0, 15, "ft")

tolerance = 5e-3
error = 1.0
#Iterate in a while loop until the MTOM error is less than 0.5% (arbitrary)
while error > tolerance:

#Get the values for m_f, m_tom, m_oem, C_D from an excel or the class II weight estimation
  CD_0_new = drag_compute(wing["mac wing"], wing["mac horizontal tail"], wing["mac vertical tail"], wing["Wing area"], wing["area horizontal tail"], wing["area vertical tail"], wing["diameter engine"], l_engine, wing["Fuselage diameter"], wing["L1"], wing["L2"], wing["L3"], wing["Mach number cruise"], rho_cruise, T_cruise, wing["dynamic viscosity cruise"])
  m_oem_new = final_weight_computation(m_tom_old, mzfw_old, wing["area horizontal tail"], wing["Wing span"], wing["Root chord"]*wing['t/c MS0313'], sweep_wing_halfchord, n_ult, wing["area horizontal tail"], wing["area vertical tail"], wing['sweep H'], wing['sweep V'], V_dive, wing['H tail arm'], wing['Fuselage diameter'], wing['Fuselage diameter'], math.pi*(wing['Fuselage diameter']**2)/4)
  mzfw_new = m_oem_new + wing['Maximum payload']
  #The m_f can be obtained with the range of the aircraft, however this may give issues in the convergence of the code since the aircraft is becoming heavier while the dimensional components are not growing
  #For that reason reitaration of the aircaft dimensional variables is probably needed
  #The fuel of the aircraft has to fit in the wings, there is a function to get the wing volume in a previous file
  m_f_new = ((math.exp(wing['Design range']*9.81*wing['Thrust specific fuel consumption']/(wing['Cruise speed']*wing['Lift-to-drag cruise'])))-1)*mzfw_new
  m_tom_new = mzfw_new + m_f_new
  error = abs(m_tom_new - m_tom_old)/m_tom_old 

#the update of the excel was deleted since more testing is probably needed before the final reliable values are achieved
  m_tom_old = m_tom_new
  m_oem_old = m_oem_new
  m_f_old = m_f_new
  mzfw_old = mzfw_new
  CD_0_old = CD_0_new
  

print(f"The final values are m_f = {m_f_new} kg\nm_tom_new = {m_tom_new} kg\nm_oem = {m_oem_new} kg\nC_D0 = {CD_0_new}")



