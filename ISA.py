import math

g0 =  9.80665 # [m/s^2]
R = 287.05287 # [J/kgK]
P_STA_SL = 101325 # [Pa]
T_STA_SL = 288.15 # [K]
rho_SL = P_STA_SL/(R*T_STA_SL) # [kg/m3]
ft = 0.3048 # [m]
FL = 100.0 * ft # [m]

layers = [[0,11000,-0.0065], # [min_h, max_h, a]
          [11000,20000,0], 
          [20000,32000,0.0010], 
          [32000,47000,0.0028], 
          [47000,51000,0], 
          [51000,71000,-0.0028], 
          [71000,86000,-0.0020]]

def find_ISA(data_array):
    [h0,h1,a,T0,P0] = data_array
    T1 = T0 + a*(h1-h0)
    if a == 0:
        P1 = P0 * math.exp(-g0/(R*T1)*(h1-h0))
    else:
        P1 = P0*(T1/T0)**(-g0/(a*R))
    rho_1 = P1/(R*T1)
    return [T1, P1, rho_1]

def calculate(h1, h0, T0, unit = "meters"):
    global T_SL
    if(h0 > 11000):
        return ValueError

    if(unit == "meters"):
        multiplier = 1
    elif(unit == "ft"):
        multiplier = ft
    elif(unit == "FL"):
        multiplier = FL

    h1 = multiplier * h1

    if(h1 > 86000):
        print("Sorry, I can only do altitudes up to 86000m")
    else:
        T0 = 273.15 + T0
        T_SL = T0 - 0.0065*(h0-0)

        # re_assign vars as start of ISA calc
        T0 = T_SL
        P0 = P_STA_SL # always 101325 Pa according to ISA

        for i in layers:
            if(i[1] < h1):
                [T0,P0,rho_0] = find_ISA([*i, T0, P0])
            else:
                [T1, P1, rho_1] = find_ISA([i[0], h1, i[2], T0, P0])
                break

        return [P1, T1, rho_1]
