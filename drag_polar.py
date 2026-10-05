import math
import numpy as np
import pandas as pd

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

def calculate_C_L_alpha(C_l_alpha: float, AR: float, M: float, Lambda_05C_rad: float, eta: float = 0.95,) -> float:
    """Calculates the 3D wing lift curve slope (C_L_alpha) based on the DATCOM formula. (not valid for airfoils)

    Parameters:
    ----------
    C_l_alpha : float
        2D airfoil lift curve slope (per radian, usually ~2 * pi).
    AR : float
        Wing aspect ratio.
    M : float
        Mach number (must be < 1.0 for subsonic flow).
    Lambda_05C_rad : float
        Sweep angle of the mid-chord line (0.5c) in RADIANS.
    eta : float, optional
        Airfoil efficiency factor (default is 0.95).

    Returns:
    -------
    float
        3D wing lift curve slope (C_L_alpha) per radian.
    """
    if M >= 1.0:
        raise ValueError("Mach number M must be subsonic (< 1.0).")

    # Compressibility correction factor
    beta = math.sqrt(1 - M**2)

    # Calculate inner sweep term
    tan_lambda = math.tan(Lambda_05C_rad)
    sweep_term = 1 + (tan_lambda**2) / (beta**2)

    # Calculate main radical term
    radical_term = 4 + ((AR * beta) / eta) ** 2 * sweep_term

    # Compute C_L_alpha
    C_L_alpha = (C_l_alpha * AR) / (2 + math.sqrt(radical_term))

    return C_L_alpha

def calculate_C_dw(M: float, M_DD: float) -> float:
    """Calculates wave drag coefficient C_d,w based on Mach number and Drag Divergence Mach number.

    Parameters:
    ----------
    M : float
        Current freestream Mach number.
    M_DD : float
        Drag divergence Mach number.

    Returns:
    -------
    float
        Wave drag coefficient (C_d,w).
    """
    # Calculate denominator inside the brackets: 1 + 2.5 * (M_DD - M) / 0.05
    # Simplifying 2.5 / 0.05 gives 50 * (M_DD - M)
    bracket_term = 1.0 + 2.5 * (M_DD - M) / 0.05

    # Prevent division by zero if bracket term is 0
    if bracket_term == 0:
        raise ZeroDivisionError(
            "Bracket term resulted in zero denominator in the exponent."
        )

    # Compute C_d,w using exponent -1
    C_dw = 0.002 * (bracket_term**-1)

    return C_dw
def calculate_alpha_c(C_L_alpha: float, alpha_CL0: float, C_L_design: float) -> float:
  """
  Parameters:
  -------------
  C_L_alpha: alpha slope of the aircraft
  alpha_CL0: angle of attack for which lift is 0
  C_L_design: CL that the aircraft have to achieve (not the same as the wing)
  
  """
  alpha_c = alpha_CL0 + C_L_alpha/C_L_design
  return alpha_c

def calculate_alpha_s(alpha_CL0: float, CL_alpha: float, delta_alpha_CLMAX:float = 1.75, CL_max:float = 1.61955 ) -> float:
    alpha_s = CL_max/CL_alpha + alpha_CL0 + delta_alpha_CLMAX
    return alpha_s

def calculate_e(AR:float, sweep:float = 0)->float:
    e1 = 2/(2-AR+math.sqrt(4+(AR^2)*(1+math.tan(sweep)^2)))
    e2 = 1.78*(1-0.045*AR^(0.68))
    return e1, e2
#COMPUTE AREA

title = "MS(1)-0313.txt"
def airfoil_unitarea(title):
  df = pd.read_csv(title, sep=r"\s+", skiprows=0, header=None)
  LE = df[0].to_numpy()
  TE = df[1].to_numpy()
  A= 0.5 * np.abs(np.dot(LE, np.roll(TE, 1)) - np.dot(TE, np.roll(LE, 1)))
  return A
def calculate_volume(title, AR:float, S:float, taperRatio:float):
  df = pd.read_csv(title, sep=r"\s+", skiprows=0, header=None)
  LE = df[0].to_numpy()
  TE = df[1].to_numpy()
  A1 = 0.5 * np.abs(np.dot(LE, np.roll(TE, 1)) - np.dot(TE, np.roll(LE, 1)))
  b = math.sqrt(b*AR)
  cr = 2*S/(b*(1+taperRatio))
  V = 0
  dy = 0.01

  for i in range(b/dy):
    c = cr + (taperRatio-1)*cr*y*2/b
    A = A1*c^2
    V
  
  
  

#CL max is just clmax*0.9 which is 1.61955
#M_DD is gotten from a supplier graph
#Cd airfoil is gotten from XFLR5
                                    