import pandas as pd
import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import ISA

# import python files as modules
import flaps
import matching_diagram
import CL_alpha
import mass_estimation
import drag_polar


# 1 - get preliminary wing planform design

df = pd.read_excel('wing_initial_planform.xlsx')
wing_unclean = df.set_index('Parameter')['value'].to_dict()
wing = {}
for key, value in wing_unclean.items():
    try:
        wing[key] = float(value)
    except (ValueError, TypeError): # for strings or blanks
        wing[key] = value  

# print(wing)

# 2 - append all design choices

wing["choices"] = {
    # "flap_type": "single slotted fowler flap"
    "flap_type": "single slotted flap with slats"
}

# 3 - airfoil data input 

# inputted from airfoil analysis
wing["Aspect ratio"] = 8
wing.update({  
    "Cl max clean": 1.7995,           
    "a0L": -3.05,   # [deg]           
})



# mass estimation
[mtom, m_oem, m_f] = mass_estimation.get_values(wing)
print(mtom, m_oem, m_f)
wing["Fuel mass"] = m_f
wing["MTOM"] = mtom
wing["OEM"] = m_oem

# matching diagram
[design_point, S, b, T] = matching_diagram.get_wing_graph(wing) # [[W_S, T_W], S, b, T]
wing["Thrust-to-weight"] = design_point[1]
wing["Wing loading (initial, before any fuel has been burned)"] = design_point[0]
wing["Total thrust"] = T*1000
wing["Wing area"] = S
wing["Wing span"] = b

CL_alpha_values = CL_alpha.calculate(wing)
wing.update(CL_alpha_values)
# print(CL_alpha_values)

[CL_alpha_flapped, a0L_L, a0L_TO] = flaps.flaps_calculation(wing)
print(CL_alpha_flapped, a0L_L, a0L_TO)
wing["CL alpha flapped"] = CL_alpha_flapped
wing["a0L_L"] = a0L_L
wing["a0L_TO"] = a0L_TO

wing_volume = drag_polar.calculate_volume("MS(1)-0313.txt", wing["Aspect ratio"], S, wing["Taper ratio"], 0.2, 0.8)
fuel_volume = m_f/805 
print(wing_volume, fuel_volume) 



# 4 - iterative loop for maximizing SAR

# assuming zero sweep, and most optimal taper to get elliptical
SAR_max = 0
best_mtom = float('inf')
best_AR = 0
optimal_wing = wing.copy()
ar_sar_arr = []

for AR in np.arange(7,13,0.1): # [-]
    wing["Aspect ratio"] = AR

    # design point creation
    [mtom, m_oem, m_f] = mass_estimation.get_values(wing)
    wing["Fuel mass"] = m_f
    wing["MTOM"] = mtom
    wing["OEM"] = m_oem

    wing_loading_difference = 0
    wing_condition = False
    while(not wing_condition):
        [design_point, S, b, T] = matching_diagram.get_wing_graph(wing, wing_loading_difference) # [[W_S, T_W], S, b, T]
        wing["Thrust-to-weight"] = design_point[1]
        wing["Wing loading (initial, before any fuel has been burned)"] = design_point[0]
        wing["Total thrust"] = T*1000
        wing["Wing area"] = S
        wing["Wing span"] = b

        t = 0.4 # [-] 
        c_r = 2*S/((1+t)*b)
        c_t = t*c_r
        sweep_LE = math.atan(1/4 * 2*c_r/b * (1-t)) / math.pi * 180
        sweep_1_2 = math.atan(math.tan(sweep_LE - 1/4 * 2*c_r/b * (1-t)))/math.pi * 180 # [deg]

        b = math.sqrt(AR * S)
        sweep_1_2 = math.atan(math.tan(sweep_LE - 1/2 * 2*c_r/b * (1-t)))/math.pi * 180 # [deg]
        y_MAC = b/6 * (1 + 2*t)/(1 + t)
        MAC = 2/3 * c_r * (1 + t + t**2)/(1 + t)
        # print(f"wing area {S}")

        wing["Sweep quarter chord"] = 0 
        wing["Taper ratio"] = 0.4
        wing["Wing span"] = b
        wing["Tip chord"] = c_t
        wing["Root chord"] = c_r
        wing["Spanwise position MAC"] = y_MAC
        wing["MAC"] = MAC
        wing["Sweep leading edge"] = sweep_LE
        wing["Oswald factor"] = drag_polar.calculate_e(AR, sweep_1_2)

        
        wing_volume = drag_polar.calculate_volume("MS(1)-0313.txt", AR, S, t, 0.3, 0.8)
        fuel_volume = m_f/805 #[m^3]   

        if(wing_volume*0.8>=fuel_volume):
            wing_condition = True

        wing_loading_difference += 50
        # print(S, wing_volume*0.8-fuel_volume)



    CL_alpha_values = CL_alpha.calculate(wing)
    # wing.update(CL_alpha_values)
    # print(CL_alpha_values)

    [CL_alpha_flapped, a0L_L, a0L_TO] = flaps.flaps_calculation(wing)
    # print(CL_alpha_flapped, a0L_L, a0L_TO)
    wing["CL alpha flapped"] = CL_alpha_flapped
    wing["a0L_L"] = a0L_L
    wing["a0L_TO"] = a0L_TO

    # print(wing_volume)

    # Calculating SAR
    C_T = 22*wing["Bypass ratio"] ** -0.19
    V_CR = wing["Cruise speed"]
    K = 1/(math.pi * AR * wing["Oswald factor"])

    [P_cruise, T_cruise, rho_cruise] = ISA.calculate(41000, 0, 15, "ft")
    W1 = wing["MTOM"]*9.80665
    W_fuel = wing["Fuel mass"]*9.80665
    W2 = W1 - W_fuel
    q = 0.5 * rho_cruise * V_CR**2
    CL_des = 1.1 * (1/q) * (W1 + W2)/(2 * wing["Wing area"])

    M_DD = 0.95 - 0.131 - CL_des/10
    Delta_CD_w = drag_polar.calculate_C_dw(wing["Mach number cruise"], M_DD)
    CD = wing["CD0"] + K*CL_des**2 + Delta_CD_w
    D = CD * 1/2 * rho_cruise * V_CR**2

    SAR = V_CR/(D*C_T)
    # if(wing_volume <= fuel_volume):
    #     pass
    ar_sar_arr.append([float(AR), mtom, S, float(SAR), D, CL_des]) 
    if(SAR >= SAR_max and mtom <= best_mtom):
        SAR_max = SAR
        best_AR = AR
        optimal_wing = wing.copy()
    if(round(AR,1) == 11):
        print(AR, SAR, S, b, mtom)
        final_wing = wing.copy()
    


# print(SAR_max, best_AR)
# ar_sar_arr = sorted(ar_sar_arr, key=lambda x: x[3], reverse=True)
# for i in range(len(ar_sar_arr)):
#     print(ar_sar_arr[i])

df = pd.DataFrame(final_wing)
df.to_excel("final.xlsx")
# print(final_wing)
