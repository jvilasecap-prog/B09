import math
import numpy as np
import ISA



def CL_cl_ratio(sweep):
    # calculated sharpness was above 2.5, hence using that respective graph
    # using linear approximation up to 27.5 degrees sweep

    x_graph = [0, 27.5] # [deg]
    y_graph = [0.9, 0.8] # [-]

    return np.interp(sweep, x_graph, y_graph)

def d_alpha_stall_calculation(sweep):
    # calculated sharpness were rounded to be above 4 at all times, hence using that respective graph
    # using quadratic approximation using approximation of points take from graph by hand

    x_graph = [0, 10, 20, 30, 40, 50, 60] # [deg]
    y_graph = [2.2, 2, 2.1, 2.3, 2.5, 2.8, 3.3] # [deg]

    # ax^2 + bx + c
    [a, b, c] = np.polyfit(x_graph, y_graph, 2)

    return (a*sweep**2 + b*sweep + c)



def calculate(wing):

    # ASSUMED


    # 1 - Cl required

    [P_cruise, T_cruise, rho_cruise] = ISA.calculate(41000, 0, 15, unit = "ft")
    print([P_cruise, T_cruise, rho_cruise])
    SoS = math.sqrt(1.4 * ISA.R * T_cruise)
    V_CR = wing["Mach number cruise"] * SoS
    
    W1 = wing["MTOM"]*9.80665
    W_fuel = wing["Fuel mass"]*9.80665
    W2 = W1 - W_fuel
    q = 0.5 * rho_cruise * V_CR**2
    CL_des = 1.1 * (1/q) * (W1 + W2)/(2 * wing["Wing area"])
    Cl_req = CL_des/(math.cos(math.radians(wing["Sweep leading edge"])))
    # print(CL_des, Cl_req)
    


    # 2 - CL_alpha using semi-empiric formula

    AR = wing["Aspect ratio"]
    beta = math.sqrt(1 - wing["Mach number cruise"])
    eta = 0.95

    #   calculation for half sweep angle
    Sweep_quarter_chord = math.radians(wing["Sweep quarter chord"])
    t = wing["Taper ratio"]
    half_sweep_angle = math.atan(math.tan(Sweep_quarter_chord - 4*(0.5 - 0.25)/AR * (1-t)/(1+t)))

    #   DATCOM CL alpha [rad^-1] and trim angle
    CL_alpha = (2 * math.pi * AR)/(2 + math.sqrt(4 + (AR*beta/eta)**2 * (1 + (math.tan(half_sweep_angle)/beta)**2)))
    atrim = CL_des/CL_alpha + wing["a0L"]*(math.pi/180) #[deg]
    # print(CL_alpha)
    


    # 3 - CL_max and stall angle

    CLmax_Clmax_ratio = CL_cl_ratio(wing["Sweep leading edge"])
    Clmax = wing["Cl max clean"]
    CL_max = CLmax_Clmax_ratio * Clmax 

    d_alpha_stall = d_alpha_stall_calculation(wing["Sweep leading edge"])
    a0L = wing["a0L"]
    alpha_stall = CL_max/CL_alpha + a0L + d_alpha_stall

    return {
        "CL max clean": CL_max,
        "CL alpha clean": CL_alpha,
        "Trim angle": atrim,
        "alpha stall": alpha_stall 
    }