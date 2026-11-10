import pandas as pd

#import python files as modules
import mass_estimation
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
tolerance = 5e-3
error = 1.0
#Iterate in a while loop until the MTOM error is less than 0.5% (arbitrary)
while error > tolerance:

#Get the values for m_f, mtom, m_oem, C_D from an excel or the class II weight estimation
  mtom_new, m_oem, m_f, C_D0 = #here classII weight estimation would be called and the old values would be used as input

  error = abs(mtom_new - mtom_old)/mtom_old
#Update dictionary with new values for next iteration
  wing['Fuel mass'] = m_f
  wing['MTOM'] = mtom_new
  wing['OEM'] = m_oem
  wing['CD_0'] = C_D0

  mtom_old = mtom_new

print(f"The final values are m_f = {wing['Fuel mass']} kg\nmtom_new = {wing['MTOM']} kg\nm_oem = {wing['OEM']} kg\nC_D0 = {C_D0}")



