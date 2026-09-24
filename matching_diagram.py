import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
import ISA


# assumed values:
takeoff_landing_altitude = 0 # [m]
takeoff_landing_temperature = 15 # [*C]
C_LFL = 0.45 # (ADSEE)
theta_t_break = 1.08 # (ADSEE)
AR = 8.5 

climb_rate_requirement = 0.5 # [m/s]
climb_rate_requirement_altitude = 12496 + 300 # [m] little above cruise altitude

clearance_height = 11 # [m] 
k_T = 0.85

# semi_empirical calculated values (change with airfoil calc)
C_D0 = 0.0212 # at cruise
e = 0.812212844 # at cruise

climb_gradient_L_D_dict = {
    "119":{
        "climb gradient": 3.2,
        "C_D0": 0.0842,
        "e": 0.903212844,
        "N_e": 2
    },
    "a": {
        "climb gradient": 0,
        "C_D0": 0.0582,
        "e": 0.851212844,
        "N_e": 1
    },
    "b": {
        "climb gradient": 2.4,
        "C_D0": 0.0407,
        "e": 0.851212844,
        "N_e": 1
    },
    "c": {
        "climb gradient": 1.2,
        "C_D0": 0.0212,
        "e": 0.812212844,
        "N_e": 1
    },
    "d": {
        "climb gradient": 2.1,
        "C_D0": 0.0667,
        "e": 0.903212844,
        "N_e": 1
    }
}

def switch_case_constraint(constraint, wing, x):
    match constraint:
        case "Take-off field length": 
            return takeoff_length(wing, x)
        case "Cruise Speed": 
            return cruise_speed(wing, x)
        case "Climb rate": 
            return climb_rate(wing, x)
        case "Climb gradient CS25.119": 
            return climb_gradient(wing, x, "119")
        case "Climb gradient CS25.121a": 
            return climb_gradient(wing, x, "a")
        case "Climb gradient CS25.121b": 
            return climb_gradient(wing, x, "b")
        case "Climb gradient CS25.121c": 
            return climb_gradient(wing, x, "c")
        case "Climb gradient CS25.121d": 
            return climb_gradient(wing, x, "d")
        case _:
            return None

def alpha_t_calc(wing, altitude, W_S, M=0, C_L=0):
    [P1, T1, rho_1] = ISA.calculate(altitude, takeoff_landing_altitude, takeoff_landing_temperature)
    if(M == 0):
        V = math.sqrt(W_S * 2/rho_1 * 1/C_L)
        a = math.sqrt(1.4 * ISA.R * T1)
        M = V/a

    # thrust laspse effect on temp and pressure
    T_t = T1*(1 + 0.2*M**2)
    P_t = P1*(1 + 0.2*M**2)**3.5

    del_t = P_t/ISA.P_STA_SL
    theta_t = T_t/ISA.T_SL

    B = wing["Bypass ratio"]

    # alpha_t dependant calculation:

    if(0 < B < 5 and theta_t < theta_t_break):
        alpha_t = del_t
    elif(0 < B < 5 and theta_t >= theta_t_break):
        alpha_t = del_t*(1 - 2.1 * (theta_t-theta_t_break)/theta_t)
    elif(5 < B < 15 and theta_t < theta_t_break):
        alpha_t = del_t*(1 - (0.43 + 0.014*B)*math.sqrt(M))
    elif(5 < B < 15 and theta_t >= theta_t_break):
        alpha_t = del_t*(1 - (0.43 + 0.014*B)*math.sqrt(M) - 3*(theta_t-theta_t_break)/(1.5+M))

    return alpha_t

def minimum_speed(wing):
    VS0 = wing["Approach speed"] / 1.23
    rho = wing["Density sea level"]
    beta = wing["Mass fraction landing"]
    CL_max_landing = wing["CL max landing"]

    return rho/(2*beta) * VS0**2 * CL_max_landing

def landing_length(wing):
    beta = wing["Mass fraction landing"]
    L_LF = wing["Landing distance"]
    rho = wing["Density sea level"]
    CL_max_landing = wing["CL max landing"]

    return (1/beta) * (1/0.6) * (L_LF/C_LFL) * (rho * CL_max_landing / 2)

def takeoff_length(wing, W_S):
    [P1, T1, rho_1] = ISA.calculate(0, takeoff_landing_altitude, takeoff_landing_temperature)

    L_D_values = climb_gradient_L_D_dict["a"]
    e = L_D_values["e"]
    C_D0 = L_D_values["C_D0"]
    C_L_req = math.sqrt(C_D0 * math.pi * AR * e)
    N_e = wing["Number of engines"]
    h_2 = clearance_height
    L_TO = wing["Take-off distance"]


    alpha_t = alpha_t_calc(wing, 0, W_S, C_L=C_L_req )

    
    T_W =  (1/alpha_t) * 1.15*(math.sqrt((N_e/(N_e-1)) * W_S/(L_TO * rho_1 * k_T * 9.80665 * math.pi * AR * e)) + (N_e/(N_e-1)) * 4*h_2/L_TO)
    return T_W # placeholder for now

def cruise_speed(wing, W_S):
    [P1, T1, rho_1] = ISA.calculate(wing["Cruise altitude IN METERS"], takeoff_landing_altitude, takeoff_landing_temperature)
    SoS = math.sqrt(1.4 * ISA.R * T1)

    V_CR = wing["Mach number cruise"] * SoS
    beta = wing["Mass fraction cruise"]
    alpha_t = alpha_t_calc(wing, wing["Cruise altitude IN METERS"], W_S, wing["Mach number cruise"])

    term_1 = (C_D0 * 1/2 * rho_1 * V_CR**2)/(beta * W_S)
    term_2 = (beta * W_S)/(math.pi * AR * e * 1/2 * rho_1 * V_CR**2)
    T_W = beta/alpha_t*(term_1 + term_2)

    return T_W 

def climb_rate(wing, W_S):
    [P1, T1, rho_1] = ISA.calculate(climb_rate_requirement_altitude, takeoff_landing_altitude, takeoff_landing_temperature)

    C_L_req = math.sqrt(C_D0 * math.pi * AR * e) # ADSEE Book (pg 162)
    beta = wing["Mass fraction cruise"]
    alpha_t = alpha_t_calc(wing, climb_rate_requirement_altitude, W_S, C_L=C_L_req )

    term_1 = math.sqrt(climb_rate_requirement**2/(beta*W_S) * rho_1/2 * math.sqrt(C_D0 * math.pi * AR * e))
    term_2 = 2*math.sqrt(C_D0/(math.pi * AR * e))
    T_W = beta/alpha_t*(term_1 + term_2)
    return T_W 

def climb_gradient(wing, W_S, type):
    [P1, T1, rho_1] = ISA.calculate(0, takeoff_landing_altitude, takeoff_landing_temperature)

    L_D_values = climb_gradient_L_D_dict[type]
    e = L_D_values["e"]
    C_D0 = L_D_values["C_D0"]
    C_L_req = math.sqrt(C_D0 * math.pi * AR * e)
    climb_grad = L_D_values["climb gradient"]/100

    N_e = L_D_values["N_e"]
    beta = wing["Mass fraction take-off"]
    alpha_t = alpha_t_calc(wing, 0, W_S, C_L=C_L_req )

    
    T_W = (wing["Number of engines"]/N_e) * (beta/alpha_t) * (climb_grad + 2*math.sqrt(C_D0/(math.pi*AR*e)))
    return T_W # placeholder for now

def plot_matching_diagram(x, y, design_point, landing_length_W_S, minimum_speed_W_S):
    plt.figure(figsize=(10, 6))

    for constraint, y_vals in y.items():
        plt.plot(x, y_vals, label=constraint, linewidth=2)

    plt.scatter(
        [design_point[0]],
        [design_point[1]],
        color="red",
        zorder=5,
        label="Design Point",
    )

    plt.axvline(x = landing_length_W_S, label = 'Minimum speed')
    plt.axvline(x = minimum_speed_W_S, label = 'Landing field length')

    plt.xlim(0, max(x))
    plt.ylim(0,design_point[1] * 1.5) # add enough space for y-scale

    plt.xlabel("Wing Loading ($W/S$) [N/m²]")
    plt.ylabel("Thrust-to-Weight Ratio ($T/W$) [-]")
    plt.title("Aircraft Matching Diagram")
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.tight_layout()
    plt.show()

def get_wing_graph(wing):

    # 1 - get vertical constraints

    minimum_speed_W_S = minimum_speed(wing)
    landing_length_W_S = landing_length(wing)

    # 2 - get remaining constraints points against w/s to plot

    x = list(range(100,round(max(minimum_speed_W_S, landing_length_W_S))+200, 20)) # adaptive list for x (don't start from 0)
    y = {
        "Take-off field length": [],
        "Cruise Speed": [],
        "Climb rate": [],
        "Climb gradient CS25.119": [],
        "Climb gradient CS25.121a": [],
        "Climb gradient CS25.121b": [],
        "Climb gradient CS25.121c": [],
        "Climb gradient CS25.121d": [],
    }

    for i in x:
        for constraint in list(y.keys()):
            y[constraint].append(switch_case_constraint(constraint, wing, i))

    # 3 - get design point

    design_point = [min(minimum_speed_W_S, landing_length_W_S), 0]
    for constraint in list(y.keys()):
        val = switch_case_constraint(constraint, wing, design_point[0])
        if val > design_point[1]:
            design_point[1] = val

    # 4 - get required calculations for wings and thrust

    S = wing["MTOM"]*9.80665/design_point[0] # [m^2]
    b = math.sqrt(wing["Aspect ratio"]*S) # [m]
    T = design_point[1]*9.80665*wing["MTOM"]/1000 # [kN]

    # 5 - plot and print everything

    plot_matching_diagram(x, y, design_point, landing_length_W_S, minimum_speed_W_S)

    output = [design_point, S, b, T]

    return output


# TESTING

# df = pd.read_excel('wing_initial_planform.xlsx')
# wing_unclean = df.set_index('Parameter')['value'].to_dict()
# wing = {}
# for key, value in wing_unclean.items():
#     try:
#         wing[key] = float(value)
#     except (ValueError, TypeError):
#         wing[key] = value  

# print(climb_gradient(wing,300,"a"))
