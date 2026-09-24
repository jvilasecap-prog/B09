"""
Module for all calculations connected to ailerons
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

EXCEL = Path(__file__).parent.resolve() / "wing_initial_planform.xlsx"

EPSILON = 1e-9

# Wing data
LIFT_CURVE_SLOPE_APP = 2 * np.pi  # per radian in landing configuration
LIFT_CURVE_SLOPE = 0.186 * 57.296  # per radian in cruise

# Aileron data
CHORD_FRACTION = 0.22
AILERON_START_FRAC = 0.6
AILERON_END_FRAC = 0.86
MAX_DEFLECTION = np.pi / 12


wing_data = pd.read_excel(EXCEL, usecols="A:B")

wing_data = wing_data.set_index(wing_data.columns[0])[wing_data.columns[1]].to_dict()

for key, value in wing_data.items():
    try:
        wing_data[key] = float(wing_data[key])
    except ValueError:
        # print(f"Value '{value}' with key '{key}' could not be converted")
        pass


def wing_shape_calc() -> tuple[float, float]:
    """
    Used to calculate the trapezoidal wing shape for integration.
    Assuming the trailing edge is equivalent to the x-axis, the leading edge
    can be modelled as a linear function to obtain a trapezoid.

    Returns a and b, coefficients for the leading edge.
    """

    tip = wing_data["Wing Span"] / 2
    if tip < EPSILON:
        raise ValueError("Tip cannot be at the same place as root.")

    a = (wing_data["Tip chord"] - wing_data["Root chord"]) / tip
    b = wing_data["Root chord"]

    return a, b


def rolling_moment_coeff_calc(a: float, b: float) -> float:
    """
    Implements the yellow formula from ADSEE lecture 3 slide 57.
    a, b - shape of the wing according to wing_shape()
    """

    def antiderivative(y: float) -> float:
        return (a / 3) * y**3 + (b / 2) * y**2

    factor = (2 * LIFT_CURVE_SLOPE * aileron_efficiency_calc(CHORD_FRACTION)) / (
        wing_data["Wing area"] * wing_data["Wing Span"]
    )
    y_1 = AILERON_START_FRAC * wing_data["Wing Span"] / 2
    y_2 = AILERON_END_FRAC * wing_data["Wing Span"] / 2
    return factor * (antiderivative(y_2) - antiderivative(y_1))


def roll_damping_coeff_calc(a: float, b: float) -> float:
    """
    Implements the yellow formula from ADSEE lecture 3 slide 59.
    a, b - shape of the wing according to wing_shape()
    """
    y_max = wing_data["Wing Span"] / 2

    factor = -(4 * (LIFT_CURVE_SLOPE + wing_data["Zero-lift drag coefficient"])) / (
        wing_data["Wing area"] * (wing_data["Wing Span"] ** 2)
    )

    return factor * (a / 4 * y_max**4 + b / 3 * y_max**3)


def aileron_efficiency_calc(
    chord_fraction: float | list | np.ndarray,
) -> float:
    """
    Calculate aileron efficiency based on the chord fraction the aileron
    occupies.

    Parameters
    ----------
    chord_fraction : float, list, or numpy.ndarray
        Input value or values in the interval [0.0, 0.7].

    Returns
    -------
    float or numpy.ndarray
        Interpolated value or values of aileron efficiency.

    Raises
    ------
    ValueError
        If any input value is outside the graph's range.
    """
    # Approximate data points read from chord fraction / efficeincy graph
    x_data = np.array([0.0, 0.033, 0.07, 0.1, 0.19, 0.28, 0.4, 0.54, 0.7])

    tau_data = np.array([0.0, 0.1, 0.2, 0.27, 0.4, 0.5, 0.6, 0.7, 0.8])

    x_array = np.asarray(chord_fraction, dtype=float)

    if np.any((x_array < x_data[0]) | (x_array > x_data[-1])):
        raise ValueError(f"x must be between {x_data[0]} and {x_data[-1]}")

    result = np.interp(x_array, x_data, tau_data)

    # Return a Python float for a scalar input
    if result.ndim == 0:
        return float(result)

    return float(result[0])


def steady_state_roll_rate(velocity) -> float:
    wing_shape = wing_shape_calc()

    roll_damping_coeff = roll_damping_coeff_calc(*wing_shape)
    rolling_moment_coeff = rolling_moment_coeff_calc(*wing_shape)

    factor = - MAX_DEFLECTION * 2 * velocity / wing_data["Wing Span"]

    return factor * rolling_moment_coeff / roll_damping_coeff


if __name__ == "__main__":
    print(steady_state_roll_rate(60))
