import math


def drag_polar(wing):
    
    # 1 - CL_alpha using semi-empiric formula

    AR = wing["Aspect ratio"]
    beta = math.sqrt(1 - wing["Mach number cruise"])
    eta = 0.95

    ## calculation for half sweep angle
    Sweep_quarter_chord = math.radians(wing["Sweep quarter chord"])
    t = wing["Taper ratio"]
    half_sweep_angle = math.atan(math.tan(Sweep_quarter_chord - 4*(0.5 - 0.25)/AR * (1-t)/(1+t)))

    ## DATCOM [rad^-1]
    C_L_alpha = (2 * math.pi * AR)/(2 + math.sqrt(4 + (AR*beta/eta)**2 * (1 + (math.tan(half_sweep_angle)/beta)**2)))


                                    

