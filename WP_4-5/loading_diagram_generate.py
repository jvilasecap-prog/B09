"""
Main file for running strength-limited design simulations
"""

from scipy import integrate as ig
from scipy import interpolate as ip
from matplotlib import pyplot as plt
from matplotlib.patches import FancyArrowPatch
from pathlib import Path
import random as rd
import numpy as np


import aerodynamic_loads as aero

# technical constants
EPSILON = 1e-9
PRECISION = 100  # number of spanwise stations at which functions are evaluated

# physical constants
G = 9.80665  # m s^-2


# aeroplane information
WING_MASS = 658  # kg
TAPER_RATIO = 0.4  # -
FUEL_VOLUME_TOTAL = 2.9477  # m^3
WING_SPAN = 17.12222083  # m
WING_AREA = 26.6518587498352  # m^2
WING_MTOM = 7990.18
MAC_LEN = 1.651865461  # m
TANK_END = 3  # m (spanwise location of the outboard fuel tank wall)
ROOT_CHORD = 2 * WING_AREA / (WING_SPAN * ( 1 + TAPER_RATIO )) # m


# flight conditions
rho = 0.2874
V_inf = 200.6
Load_factor = 1.5

# additional information
FUEL_DENSITY = 840  # kg m^-3


def chord_dist(
    z_over_L: np.typing.ArrayLike,
    t: float = TAPER_RATIO,
) -> np.typing.NDArray[np.float64]:
    z_over_L = np.asarray(z_over_L, dtype=np.float64)

    return 1.0 - z_over_L * (1.0 - t)


def wing_weight_distribution(
    z: np.typing.ArrayLike,
    L: float = WING_SPAN / 2,
    t: float = TAPER_RATIO,
    m: float = WING_MASS / 2,
) -> np.typing.NDArray[np.float64]:
    if L <= 0:
        raise ValueError("Wing half-span must be positive.")
    if not 0 < t < 1:
        raise ValueError("Taper ratio must be between 0 and 1.")
    if m < 0:
        raise ValueError("Wing mass must be non-negative.")

    z = np.asarray(z, dtype=np.float64)

    scaling_factor = 1.0 / (L * (t + (1.0 / 3.0) * (1.0 - t) ** 2))

    distribution = scaling_factor * m * G * chord_dist(z / L, t) ** 2

    return np.where(
        (z >= 0.0) & (z <= L),
        distribution,
        0.0,
    )


def fuel_weight_distribution(
    z: np.typing.ArrayLike,
    L: float = WING_SPAN / 2,
    tank_end: float = TANK_END,
    t: float = TAPER_RATIO,
    V: float = FUEL_VOLUME_TOTAL,
    d: float = FUEL_DENSITY,
) -> np.typing.NDArray[np.float64]:
    if L <= 0:
        raise ValueError("Wing half-span must be positive.")
    if not 0 < tank_end <= L:
        raise ValueError("tank_end must lie within the interval (0, L].")
    if not 0 <= t <= 1:
        raise ValueError("t must lie within the interval [0, 1].")
    if V < 0:
        raise ValueError("V must be non-negative.")
    if d < 0:
        raise ValueError("d must be non-negative.")

    z = np.asarray(z, dtype=np.float64)

    fuel_weight_in_one_wing = 0.5 * V * d * G

    scaling_factor = 1.0 / (
        tank_end
        - tank_end**2 * (1.0 - t) / L
        + tank_end**3 * (1.0 - t) ** 2 / (3.0 * L**2)
    )

    distribution = scaling_factor * fuel_weight_in_one_wing * chord_dist(z / L, t) ** 2

    return np.where(
        (z >= 0.0) & (z <= tank_end),
        distribution,
        0.0,
    )

def get_torsion_from_lift(lift: np.typing.NDArray, ) -> np.typing.NDArray:
    spanwise_stations = lift[:,0]
    chord_distr = chord_dist(spanwise_stations / spanwise_stations[-1]) * ROOT_CHORD
    return np.column_stack((spanwise_stations, lift[:,1] * chord_distr / 4))


def internal_any(
    half_span: float = WING_SPAN / 2,
    precision: int = PRECISION,
    **kwargs: np.typing.NDArray | function | tuple[float, float],
) -> np.typing.NDArray:

    # generate the initial internal loading
    loading = np.column_stack(
        (
            np.linspace(0, half_span, precision),
            np.zeros(precision),
        )
    )

    # handle the loadings provided
    for key, val in kwargs.items():

        # determine the type of the argument considered
        if (
            isinstance(val, tuple)
            and len(val) == 2
            and isinstance(val[0], float)
            and isinstance(val[1], float)
        ):
            process_type = "point_force"
        elif isinstance(val, np.ndarray):
            distr = ip.interp1d(
                val[:, 0], val[:, 1], bounds_error=False, fill_value=0.0
            )
            process_type = "distr"
        elif callable(val):
            # this part of code is untested so try not to pass functions as an argument
            distr = val
            process_type = "distr"
        else:
            process_type = "unknown"

        match process_type:
            case "point_force":
                application_z, magnitude = val
                last_index = np.flatnonzero(loading[:, 0] <= application_z)[-1]
                loading[:last_index, 1] += magnitude
            case "distr":
                for i in range(precision):
                    loading[i, 1] += ig.quad(distr, loading[i, 0], half_span)[0]
            case _:
                raise ValueError(
                    f"The type of loading '{key}' could not be determined."
                )
    return loading


if __name__ == "__main__":

    fig, (shear, moment, loadax) = plt.subplots(1, 3)

    # Trial loadings taken from the Statics coursebook assignment 7-78
    # force_a = (0., -10.)
    # force_b = (6., -25.)
    # force_2m = (2., 15.)
    # moment_3m = (3., -20.)
    # force_distr = np.array([[4., 10.], [6., 10.]])

    # define loadings

    spanwise_stations = np.linspace(0, WING_SPAN / 2, PRECISION)

    fuel_weight = np.column_stack(
        (spanwise_stations, fuel_weight_distribution(spanwise_stations))
    )
    wing_weight = np.column_stack(
        (spanwise_stations, wing_weight_distribution(spanwise_stations))
    )
    lift_distr = -np.column_stack((-spanwise_stations, aero.lift(spanwise_stations)))

    # calculate internal forces
    shear_loading = internal_any(
        fuel_weight=fuel_weight,
        wing_weight=wing_weight,
        lift=lift_distr,
    )
    moment_loading = internal_any(
        moment_shear=np.column_stack(
            (shear_loading[:, 0], -shear_loading[:, 1])
        )  # this argument relates the shear loading to the moment, should not be changed
    )

    # draw internal forces
    shear.plot(shear_loading[:, 0], shear_loading[:, 1])
    moment.plot(moment_loading[:, 0], moment_loading[:, 1])
    loadax.plot(lift_distr[:, 0], lift_distr[:, 1])

    # format and show the graphs
    for ax in (shear, moment, loadax):
        ax.set_xlabel("z")
        ax.spines["bottom"].set_position("zero")
        ax.spines["bottom"].set_position(("data", 0))
        ax.xaxis.set_ticks_position("bottom")
        ax.xaxis.set_label_position("bottom")

    shear.set_ylabel("internal shear force")
    moment.set_ylabel("interal bending moment")
    plt.show()
