import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import ISA

# import python files as modules
import flaps
import matching_diagram
import CL_alpha
import mass_estimation


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

# 4 - iterative loop
CL_alpha_values = CL_alpha.calculate(wing)
wing.update(CL_alpha_values)
print(CL_alpha_values)

[CL_alpha_flapped, a0L_L, a0L_TO] = flaps.flaps_calculation(wing)
print(CL_alpha_flapped, a0L_L, a0L_TO)
wing["CL alpha flapped"] = CL_alpha_flapped
wing["a0L_L"] = a0L_L
wing["a0L_TO"] = a0L_TO