import math
import ISA

def get_values(wing):

    # ASSUMED
    e_p = 43*10**6 # [J/kg]
    bypass_ratio = 5.5 # (ADSEE)
    psi = wing["Psi"] # (ADSEE)
    phi = wing["Phi"] # (ADSEE)
    AR = wing["Aspect ratio"] # design choice

    C_d0 = wing["CD0"] # change based on lift curve, but first assumed with adsee book
    e = 1/(math.pi*AR*psi+ 1/phi)
    TSFC = 22*bypass_ratio**(-0.19)

    [P_cruise, T_cruise, rho_cruise] = ISA.calculate(41000, 0, 15, unit = "ft")
    print([P_cruise, T_cruise, rho_cruise])
    SoS = math.sqrt(1.4 * ISA.R * T_cruise)
    V_CR = wing["Mach number cruise"] * SoS

    jet_efficiency = (V_CR/(TSFC*10**-6))/e_p
    L_D_max_cruise = 0.5* math.sqrt((math.pi*AR*e)/(C_d0))

    R_eq_res = 805165.68 # [m] (constant adsee calc based on fixed TLARs)
    R_lost = 1/0.7*L_D_max_cruise*(12496.8+(V_CR**2/(2*9.80665)))
    R_eq = (6100*1000 + R_lost)*1.05 + R_eq_res

    m_f_mtom_ratio = 1 - math.exp((-R_eq*9.80665)/(jet_efficiency*L_D_max_cruise*e_p))
    m_oem_mtow_ratio = 0.61
    m_pl = wing["Maximum payload"]

    mtom = m_pl/(1 - m_f_mtom_ratio - m_oem_mtow_ratio)
    m_f = mtom * m_f_mtom_ratio
    m_oem = mtom * m_oem_mtow_ratio

    return [mtom, m_oem, m_f]