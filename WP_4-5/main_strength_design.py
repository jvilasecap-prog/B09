"""
Main file for running strength-limited design simulations
"""

from scipy import integrate as ig
from scipy import interpolate as ip
import random as rd
import numpy as np



# technical constants
EPSILON = 1e-9

# physical constants
G = 9.80665 # m s^-2

# aeroplane information
WING_WEIGHT = 658 # kg
TAPER_RATIO = 0.4 # -

def wing_weight_distribution(x: float, L: float, t: float = TAPER_RATIO, m: float = WING_WEIGHT) -> float:
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

    if not (0 < x < L):
        raise ValueError("Spanwise position must be between 0 and the half-span")
    if not (0 < t < 1):
        raise ValueError("Taper ratio must be between 0 and 1")

    
    def chord_dist (x):
        return 1 - (x / L) * (1 - t)

    scaling_factor = 1 / ig.quad(lambda x : chord_dist(x) ** 2, 0, 1)[0]

    return scaling_factor * m * chord_dist(x) ** 2

if __name__ == "__main__":
    for i in range (1, 10_000):
        rd.seed(42)
        l = rd.random() * 15
        w = 42000
        assert abs(ig.quad(lambda x : wing_weight_distribution(x, l, rd.random(), w), 0, l)[0] - w) < EPSILON

    