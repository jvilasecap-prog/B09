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
