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
    cf_c_ratio = 0.35/0.4 # [-] (ADSEE)
    hinge_position_fraction = 0.7 # c_hinge/c (NASA)
    Delta_CL_max_landing = 2.65 - wing["CL max clean"] # larger than CL
    Delta_CL_max_takeoff = 0.8*Delta_CL_max_landing # 0.8*CL_max_landing (ADSEE)
    Delta_a0l_airfoil_landing = -15 # [deg] (ADSEE)
    Delta_a0l_airfoil_takeoff = -10 # [deg]


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

    """1 - Delta_CL_max"""

    # calculating hinge line sweep
    Sweep_quarter_chord = math.radians(wing["Sweep quarter chord"])
    A = wing["Aspect ratio"]
    t = wing["Taper ratio"]
    Lambda_hinge_line = math.atan(math.tan(Sweep_quarter_chord - 4*(hinge_position_fraction - 0.25)/A * (1-t)/(1+t))) # (NASA)

    Swf_S_ratio = Delta_CL_max_landing/(0.9 * Delta_Cl_max * math.cos( Lambda_hinge_line ))
    # print(Swf_S_ratio)

    """2 - Delta_a0L calculation"""

    Delta_a0L_landing = Delta_a0l_airfoil_landing * Swf_S_ratio *  math.cos( Lambda_hinge_line ) # [degrees]
    Delta_a0L_takeoff = Delta_a0l_airfoil_takeoff * Swf_S_ratio *  math.cos( Lambda_hinge_line ) # [degrees]

    """3 - CL alpha flapped calculation"""
   
    S2_S1_ratio = 1 + Swf_S_ratio * (c2_c1_ratio - 1) # (ADSEE ppt)
    CL_alpha_flapped = S2_S1_ratio * wing["CL alpha clean"]

    """4 - Add contributions and create function"""

    # print(Delta_CL_max, Delta_a0L, CL_alpha_flapped)
    
    CL_max_landing = wing["CL max clean"] + Delta_CL_max_landing
    CL_max_takeoff = wing["CL max clean"] + Delta_CL_max_takeoff
    a0L_L = wing["a0L"] + Delta_a0L_landing
    a0L_TO = wing["a0L"] + Delta_a0L_takeoff

    # y = Ax + B
    A_L = CL_alpha_flapped
    B_L = -CL_alpha_flapped*a0L_L

    A_TO = CL_alpha_flapped
    B_TO = -CL_alpha_flapped*a0L_TO

    """5 - angle for CL_max_landing=2.5"""

    # print((2.5 - B_L)/A_L)



