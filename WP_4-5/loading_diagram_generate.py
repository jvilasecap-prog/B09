"""
Main file for running strength-limited design simulations
"""

from scipy import integrate as ig
from scipy import interpolate as ip
from matplotlib import pyplot as plt
import random as rd
import numpy as np



# technical constants
EPSILON = 1e-9

# physical constants
G = 9.80665 # m s^-2

# aeroplane information
WING_MASS = 658 # kg
TAPER_RATIO = 0.4 # -
FUEL_VOLUME = 2.9477 # m^3

# additional information
FUEL_DENSITY = 840 # kg m^-3


def chord_dist (x_over_L: float, t: float) -> float:
    """
    Return the local chord length normalized by the root chord length.

    The chord distribution is assumed to vary linearly along the span
    according to the taper ratio.

    Parameters
    ----------
    x_over_L : float
        Spanwise position normalized by the half-span, defined as
        :math:`x/L`, where 0 corresponds to the wing root and 1 to the
        wing tip.
    t : float
        Wing taper ratio, defined as the tip chord divided by the root
        chord (:math:`c_t/c_r`).

    Returns
    -------
    float
        Ratio of the local chord length to the root chord length,
        :math:`c/c_r`, at the specified spanwise location.

    Examples
    --------
    >>> chord_dist(0.0, 0.4)
    1.0

    >>> chord_dist(1.0, 0.4)
    0.4

    >>> chord_dist(0.5, 0.4)
    0.7
    """
    return 1 - (x_over_L) * (1 - t)

def wing_weight_distribution(x: float, L: float, t: float = TAPER_RATIO, m: float = WING_MASS) -> float:
    """
    Evaluate the mass distribution function of a half-wing at spanwise position `z`,
    assuming constant density.

    The half-wing has half-span `L` and taper ratio `t`.

    Parameters
    ----------
    z : float
        Spanwise position along the half-wing, measured from the wing root.
        Must lie within the inclusive interval [0, L].
    L : float
        Half-span of the wing.
    t : float
        Wing taper ratio. Must lie within the inclusive interval [0, 1].

    Returns
    -------
    float
        Value of the half-wing mass distribution at position `z`.

    Raises
    ------
    ValueError
        If `z` is outside the interval [0, L].
        If `t` is outside the interval [0, 1].
    """

    if L <= 0:
        raise ValueError("Wing half-span must be positive")
    if not (0 < x < L):
        raise ValueError("Spanwise position must be between 0 and the half-span")
    if not (0 < t < 1):
        raise ValueError("Taper ratio must be between 0 and 1")

    # scaling factor is the area under the graph of the square of the chord distribution
    scaling_factor = 1 / (L * (t + (1 / 3) * (1 - t) ** 2))


    return scaling_factor * m * G * chord_dist(x / L, t) ** 2

def fuel_weight_distribution(x: float, L: float, tank_end: float, t: float = TAPER_RATIO, V: float = FUEL_VOLUME, d: float = FUEL_DENSITY) -> float:
    """
    Evaluate the fuel weight distribution function of a half-wing at spanwise position `x`.

    The fuel tank extends from the wing root to the spanwise position `tank_end`.
    The distribution is based on the half-span `L`, taper ratio `t`, fuel-tank
    volume `V`, and fuel density `d`. The fuel weight distribution is zero
    outside the fuel tank. The wing is trapezoidal

    Parameters
    ----------
    x : float
        Spanwise position along the half-wing, measured from the wing root.
        Must lie within the inclusive interval [0, L].
    L : float
        Half-span of the wing.
    tank_end : float
        Spanwise position of the outboard wall of the fuel tank, measured from
        the wing root. Must lie within the interval (0, L].
    t : float, optional
        Wing taper ratio. Must lie within the inclusive interval [0, 1].
        The default value is `TAPER_RATIO`.
    V : float, optional
        Fuel-tank volume. Must be non-negative. The default value is
        `FUEL_VOLUME`.
    d : float, optional
        Fuel density. Must be non-negative. The default value is
        `FUEL_DENSITY`.

    Returns
    -------
    float
        Fuel weight per unit span at position `x`. Returns zero when `x` is
        outboard of `tank_end`.

    Raises
    ------
    ValueError
        If `x` is outside the interval [0, L].
        If `tank_end` is outside the interval [0, L].
        If `t` is outside the interval [0, 1].
        If `V` is negative.
        If `d` is negative.
    """
    
    if L <= 0:
        raise ValueError("Wing half-span must be positive")

    if not 0 <= x <= L:
        raise ValueError("x must lie within the interval [0, L].")

    if not 0 < tank_end <= L:
        raise ValueError("tank_end must lie within the interval (0, L].")

    if not 0 <= t <= 1:
        raise ValueError("t must lie within the interval [0, 1].")

    if V < 0:
        raise ValueError("V must be non-negative.")

    if d < 0:
        raise ValueError("d must be non-negative.")

    total_fuel_weight = V * d * G

    scaling_factor = 1 / (tank_end * (t + (1 / 3) * (1 - t) ** 2))

    return scaling_factor * total_fuel_weight * chord_dist(x / L, t) ** 2

    

if __name__ == "__main__":
    # plt.plot(np.linspace(0,1,30), chord_dist(np.linspace(0,1,30), .2))
    # plt.xlim(0,1)
    # plt.ylim(0,1)
    # plt.show()
    rd.seed(42)
    for i in range (10_000):
        l = rd.random() * 100
        te = rd.random() * l
        x = rd.random() * te
        weight = FUEL_DENSITY * FUEL_VOLUME * G
        weight_ig = ig.quad(lambda x: fuel_weight_distribution(x = x, L = l, tank_end = te), 0, te)[0]
        if not (abs(weight - weight_ig) < EPSILON):
            print (i)
            print (weight, weight_ig)
            break
        # weight_ig = ig.quad(lambda x: wing_weight_distribution(x, l), 0, l)[0]
        # assert abs(WING_MASS * G - weight_ig) < EPSILON
    print("complete")