"""
Class 2 weight estimation method by Torenbeek, for SI system
"""

import numpy as np
import pandas as pd
from pathlib import Path

EXCEL = Path(__file__).parent.resolve() / "final_wing_planform.xlsx"

EPSILON = 1e-9

ac_data = pd.read_excel(EXCEL, usecols="A:B")

ac_data = ac_data.set_index(ac_data.columns[0])[ac_data.columns[1]].to_dict()

for key, value in ac_data.items():
    try:
        ac_data[key] = float(ac_data[key])
    except ValueError:
        # print(f"Value '{value}' with key '{key}' could not be converted")
        pass


def w_struct(
    MTOW: float,
    fus_d: float,
    fus_h: float,
    fus_l: float,
    n_ult: float,
) -> float:
    return MTOW * (0.447 * np.sqrt(n_ult) * ((fus_d * fus_h * fus_l) / MTOW) ** 0.24)


def w_wing(
    MZFW: float,
    wing_area: float,
    wing_span: float,
    root_thickness: float,
    half_chord_sweep: float,
    n_ult: float,
) -> float:

    structural_span = wing_span / np.cos(half_chord_sweep)

    wing_mass = (
        6.67e-3
        * structural_span**0.75
        * (1.0 + np.sqrt(1.905 / structural_span))
        * n_ult**0.55
        * (structural_span / root_thickness)
        * (MZFW / wing_area) ** 0.30
    )

    return wing_mass


import numpy as np


def w_tailplane(
    horizontal_area: float,
    vertical_area: float,
    horizontal_sweep: float,
    vertical_sweep: float,
    dive_speed: float,
) -> float:
    """
    Calculate the combined horizontal and vertical stabilizer mass using
    Torenbeek equations (8-14) and (8-15).

    Parameters
    ----------
    horizontal_area : float
        Horizontal stabilizer area [m^2].
    vertical_area : float
        Vertical stabilizer area [m^2].
    horizontal_sweep : float
        Horizontal stabilizer sweep angle [rad].
    vertical_sweep : float
        Vertical stabilizer sweep angle [rad].
    dive_speed : float
        Design dive speed as equivalent airspeed [m/s].

    Returns
    -------
    float
        Combined stabilizer mass [kg].

    Raises
    ------
    ValueError
        If an area or dive speed is non-positive, cosine of a sweep angle
        is non-positive, or an input to f is outside [0.2, 0.45].
    """
    if horizontal_area <= 0.0 or vertical_area <= 0.0:
        raise ValueError("Stabilizer areas must be positive.")

    if dive_speed <= 0.0:
        raise ValueError("Design dive speed must be positive.")

    if np.cos(horizontal_sweep) <= 0.0:
        raise ValueError("Horizontal stabilizer sweep gives cos(sweep) <= 0.")

    if np.cos(vertical_sweep) <= 0.0:
        raise ValueError("Vertical stabilizer sweep gives cos(sweep) <= 0.")

    def f(x: float) -> float:
        if x < 0.2 or x > 0.45:
            raise ValueError(
                f"Input to f must be within [0.2, 0.45]; received {x:.4f}."
            )

        # Linear interpolation between (0.2, 10) and (0.45, 25).
        return 10.0 + (25.0 - 10.0) / (0.45 - 0.2) * (x - 0.2)

    horizontal_f_input = (
        horizontal_area**0.2
        * dive_speed
        / (1000.0 * np.sqrt(np.cos(horizontal_sweep)))
    )

    vertical_f_input = (
        vertical_area**0.2
        * dive_speed
        / (1000.0 * np.sqrt(np.cos(vertical_sweep)))
    )

    # Fixed-incidence horizontal stabilizer.
    k_h = 1.0

    # Horizontal stabilizer mounted on the fin.
    k_v = 1.0 + 0.15 * horizontal_area / vertical_area

    horizontal_mass = (
        k_h
        * f(horizontal_f_input)
        * horizontal_area
    )

    vertical_mass = (
        k_v
        * f(vertical_f_input)
        * vertical_area
    )

    return horizontal_mass + vertical_mass

def w_fuselage(
    dive_speed: float,
    tail_arm: float,
    fuselage_width: float,
    fuselage_height: float,
    gross_shell_area: float,
) -> float:
    """
    Estimate fuselage mass using Torenbeek equation (8-16), including
    corrections for pressurization and rear-mounted engines.

    All inputs use metric units, and the result is in kilograms.
    """
    if dive_speed <= 0.0:
        raise ValueError("Design dive speed must be positive.")

    if tail_arm <= 0.0:
        raise ValueError("Tail arm must be positive.")

    if fuselage_width <= 0.0 or fuselage_height <= 0.0:
        raise ValueError("Fuselage width and height must be positive.")

    if gross_shell_area <= 0.0:
        raise ValueError("Gross shell area must be positive.")

    k_wf = 0.23

    fuselage_mass = (
        k_wf
        * np.sqrt(
            dive_speed
            * tail_arm
            / (fuselage_width + fuselage_height)
        )
        * gross_shell_area**1.2
    )

    # Corrections for pressurization and rear-mounted engines.
    fuselage_mass *= 1.12

    return fuselage_mass

def w_control_surfaces(
        MTOW: float,
        k_sc: float = .64 
) -> float:
    factor = .768 #metric
    return k_sc * factor * MTOW ** (2/3)

def w_landing_gear(MTOW: float) -> float:
    """
    Estimate the combined landing-gear mass using Torenbeek equation (8-17)
    for a low-wing, jet-propelled trainer or executive aircraft.

    Parameters
    ----------
    MTOW : float
        Maximum take-off mass [kg].

    Returns
    -------
    float
        Combined main and nose landing-gear mass [kg].
    """
    if MTOW <= 0.0:
        raise ValueError("MTOW must be positive.")

    main_gear_mass = (
        15.0
        + 0.033 * MTOW**0.75
        + 0.021 * MTOW
    )

    nose_gear_mass = (
        5.4
        + 0.049 * MTOW**0.75
    )

    return main_gear_mass + nose_gear_mass


def w_propulsion_group(
    engine_mass: float = 524.4,
    number_of_engines: int = 2,
) -> float:
    """
    Estimate propulsion-group mass using Torenbeek equation (8-27)
    for pod-mounted jet engines with thrust reversers.

    Parameters
    ----------
    engine_mass : float
        Definitive dry mass of one engine [kg].
    number_of_engines : int
        Number of engines [-].

    Returns
    -------
    float
        Propulsion-group mass [kg].
    """
    if engine_mass <= 0.0:
        raise ValueError("Engine mass must be positive.")

    if number_of_engines <= 0:
        raise ValueError("Number of engines must be positive.")

    k_pg = 1.15
    k_thr = 1.18

    return k_pg * k_thr * number_of_engines * engine_mass

