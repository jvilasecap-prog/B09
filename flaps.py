import math


def c2_c1_ratio_calculation(cf_c_ratio, flap_type, c1):

    # ASSUMED (max flap deflection for landing)
    flap_deflection = 40 # [deg]

    if flap_type == "single slotted flap":
        Dc_cf_ratio = 0.15/40 * flap_deflection # ref point (0.15, 40)
    elif flap_type == "single slotted fowler flap":
        Dc_cf_ratio = 0.4 + 0.3/33 * flap_deflection # ref point (0.7,33)
        
    Dc = Dc_cf_ratio * cf_c_ratio * c1
    c2 = c1 + Dc

    return c2/c1

def flaps_calculation(wing):

    # ASSUMED VALUES:
    cf_c_ratio = 0.35/0.4 # [-]
    Swf_S_ratio = 0.3 # (TBD)
    Delta_CL_max = 2.5 - wing["CL max clean"]
    Delta_a0l_airfoil_landing = -15 # [deg]
    Delta_a0l_airfoil_takeoff = -10 # [deg]
    Lambda_hinge_line = 6/180 *math.pi # Needs calculations based on sweep angle etc.

    c2_c1_ratio = c2_c1_ratio_calculation(cf_c_ratio, "single slotted fowler flap", wing["MAC"])

    # Dictionary for values corresponding to type of flap
    flaps_dict = {
        "single slotted flap": {
            "Delta_Cl_max": 1.3,
            "max_flap_deflection": 40  # [deg]
        },
        "single slotted fowler flap": {
            "Delta_Cl_max": 1.3 * c2_c1_ratio, 
            "max_flap_deflection": 40  # [deg]
        }
    }

    flap = flaps_dict[wing["choices"]["flap_type"]]

    Delta_Cl_max = flap["Delta_Cl_max"]

    # 1 - Delta_CL_max

    Swf_S_ratio = Delta_CL_max/(0.9 * Delta_Cl_max * math.cos( Lambda_hinge_line ))
    print(Swf_S_ratio)

    # 2 - Delta_a0L

    Delta_a0L = Delta_a0l_airfoil_landing * Swf_S_ratio *  math.cos( Lambda_hinge_line )

    # 3 - CL_alpha_flapped
   
    S2_S1_ratio = 1 + Swf_S_ratio * (c2_c1_ratio - 1)
    CL_alpha_flapped = S2_S1_ratio * wing["CL_alpha_clean"]

    # 4 - Add contributions
    
    CL_max = wing["CL max clean"] + Delta_CL_max
    a0L = wing["a0L"] + Delta_a0L

    # 5 - create function CL

    # y = Ax + B
    A = CL_alpha_flapped
    B = -CL_alpha_flapped*a0L


