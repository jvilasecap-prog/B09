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

# flight conditions
rho = 0.2874
V_inf = 200.6
Load_factor = 1.5

# additional information
FUEL_DENSITY = 840  # kg m^-3


def chord_dist(z_over_L: float, t: float) -> float:
    """
    Return the local chord length normalized by the root chord length.

    The chord distribution is assumed to vary linearly along the span
    according to the taper ratio.

    Parameters
    ----------
    z_over_L : float
        Spanwise position normalized by the half-span, defined as
        :math:`z/L`, where 0 corresponds to the wing root and 1 to the
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

    return 1 - (z_over_L) * (1 - t)


def wing_weight_distribution(
    z: float, L: float, t: float = TAPER_RATIO, m: float = WING_MASS
) -> float:
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
    if not (0 <= z <= L):
        return 0.
    if not (0 < t < 1):
        raise ValueError("Taper ratio must be between 0 and 1")

    # scaling factor is the area under the graph of the square of the chord distribution
    scaling_factor = 1 / (L * (t + (1 / 3) * (1 - t) ** 2))

    return scaling_factor * m * G * chord_dist(z / L, t) ** 2


def fuel_weight_distribution(
    z: float,
    L: float,
    tank_end: float,
    t: float = TAPER_RATIO,
    V: float = FUEL_VOLUME_TOTAL,
    d: float = FUEL_DENSITY,
) -> float:
    """
    Evaluate the fuel weight distribution function of a half-wing at spanwise position `x`.

    The fuel tank extends from the wing root to the spanwise position `tank_end`.
    The distribution is based on the half-span `L`, taper ratio `t`, fuel-tank
    volume `V`, and fuel density `d`. The fuel weight distribution is zero
    outside the fuel tank. The wing is trapezoidal

    Parameters
    ----------
    z : float
        Spanwise position along the half-wing, measured from the wing root.
        If outside of [0, L], returns 0.
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
        Fuel weight per unit span at position `z`. Returns zero when `z` is
        outboard of `tank_end` or inboard of root.

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
    if not 0 <= z <= tank_end:
        return 0.
    if not 0 < tank_end <= L:
        raise ValueError("tank_end must lie within the interval (0, L].")
    if not 0 <= t <= 1:
        raise ValueError("t must lie within the interval [0, 1].")
    if V < 0:
        raise ValueError("V must be non-negative.")
    if d < 0:
        raise ValueError("d must be non-negative.")

    fuel_weight_in_one_wing = 0.5 * V * d * G

    scaling_factor = 1 / (
        tank_end
        - ((tank_end**2) * (1 - t) / L)
        + ((tank_end**3) * (1 - t) ** 2 / (3 * L**2))
    )

    return scaling_factor * fuel_weight_in_one_wing * chord_dist(z / L, t) ** 2


def draw_force(
    ax: plt.Axes,
    x: float,
    y_start: float,
    y_end: float,
    label: str | None = None,
    color: str = "tab:blue",
) -> None:
    """Draw a vertical force arrow on axis ax."""
    ax.annotate(
        "",
        xy=(x, y_end),  # Arrow tip
        xytext=(x, y_start),  # Arrow tail
        arrowprops={
            "arrowstyle": "-|>",
            "color": color,
            "linewidth": 2,
            "mutation_scale": 15,
        },
    )

    if label:
        ax.text(
            x + 0.15,
            (y_start + y_end) / 2,
            label,
            color=color,
            fontsize=12,
        )


def draw_moment(
    ax: plt.Axes,
    x: float,
    y: float,
    radius: float,
    label: str | None = None,
    clockwise: bool = False,
    color="tab:red",
) -> None:
    """Draw a curved arrow representing a moment on axes ax."""
    if clockwise:
        start = (x - radius, y)
        end = (x + radius, y)
        curvature = 1.4
    else:
        start = (x + radius, y)
        end = (x - radius, y)
        curvature = 1.4

    arrow = FancyArrowPatch(
        start,
        end,
        connectionstyle=f"arc3,rad={curvature}",
        arrowstyle="-|>",
        mutation_scale=15,
        linewidth=2,
        color=color,
    )

    ax.add_patch(arrow)

    if label:
        ax.text(
            x,
            y + 1.4 * radius,
            label,
            ha="center",
            va="bottom",
            color=color,
            fontsize=12,
        )


def read_coeff_data(file_path: str | Path) -> np.typing.NDArray[np.float64]:
    """
    Read y-span, Cl, and CmGeom data from an XFLR5 results file.

    Parameters
    ----------
    file_path : str or pathlib.Path
        Path to the XFLR5 text output file.

    Returns
    -------
    numpy.ndarray
        Array with shape (3, n), where:
        - tab[0, n] is the nth y-span value
        - tab[1, n] is the nth Cl value
        - tab[2, n] is the nth CmGeom value

    Raises
    ------
    ValueError
        If the spanwise data table cannot be found or contains no data.
    """
    rows = []
    reading_table = False

    with open(file_path, "r", encoding="utf-8") as file:
        for line in file:
            if "y-span" in line and "CmGeom" in line:
                reading_table = True
                continue

            if reading_table and "Main Wing Cp Coefficients" in line:
                break

            if reading_table:
                values = line.split()

                if len(values) >= 7:
                    try:
                        y_span = float(values[0])
                        cl = float(values[3])
                        cm_geom = float(values[6])
                    except ValueError:
                        print(
                            "every minute in the project, an exception is handled without any output to console"
                        )
                        continue

                    rows.append((y_span, cl, cm_geom))

    if not rows:
        raise ValueError("No spanwise data table was found in the file.")

    return np.asarray(rows, dtype=float)


def get_loadings_from_coeff_halfspan(
    coeff: np.typing.NDArray, rho: float, vel: float, ref_S=WING_AREA, ref_c=MAC_LEN
) -> np.typing.NDArray:
    coeff[:, 1:] *= -(1 / 2) * rho * vel**2 * ref_S
    coeff[:, 2] *= ref_c
    return coeff[coeff[:, 0] > 0.0]


def shear_force_fbd(
    ax: plt.Axes,
    half_span: float = WING_SPAN / 2,
    **kwargs: np.typing.NDArray | tuple[function, float, float] | tuple[float, float],
) -> None:
    """ """
    # Plot the wing as a cantilevered beam

    ax.plot(
        [0, half_span],
        [0, 0],
        color="black",
        linewidth=10,
        solid_capstyle="butt",
    )
    ax.plot(
        [0, 0],
        [-half_span / 20, half_span / 20],
        color="black",
        linewidth=6,
        solid_capstyle="butt",
    )

    # info for limits of the axis
    largest_loading_value = 0

    # handle the loadings provided
    for key, val in kwargs.items():

        # determine the type of the argument considered
        if isinstance(val, tuple):
            if len(val) == 3 and callable(val[0]):
                process_type = "distr_func"
            elif len(val) == 2 and isinstance(val[0], float):
                process_type = "point_force"
            else:
                process_type = "unknown"

        elif isinstance(val, np.ndarray):
            process_type = "distr_precomp"

        else:
            process_type = "unknown"

        # plot the forces and distributed forces
        match process_type:
            case "point_force":
                draw_force(
                    ax=ax,
                    x=val[0],
                    y_start=0,
                    y_end=val[1],
                    label=key,
                )
                largest_loading_value = max(largest_loading_value, abs(val[1] * 1.1))
            case "distr_precomp":
                n_rows = val.shape[0]
                for i in np.linspace(0, n_rows - 1, 10):
                    e = round(i)
                    draw_force(ax=ax, x=val[e, 0], y_start=0, y_end=val[e, 1])
                ax.plot(val[:, 0], val[:, 1])
                largest_loading_value = max(
                    largest_loading_value, np.max(np.abs(val[:, 1]))
                )
            case "distr_func":
                data_z = np.linspace(val[1], val[2], 100)
                data_y = np.array([val[0](z) for z in data_z])
                for i in range(0, 100, 10):
                    draw_force(ax=ax, x=data_z[i], y_start=0, y_end=data_y[i])
                ax.plot(data_z, data_y)
                largest_loading_value = max(
                    largest_loading_value, np.max(np.abs(data_y))
                )
            case _:
                raise ValueError(
                    f"The type of loading '{key}' could not be determined."
                )

        ax.set_ylim(-largest_loading_value, largest_loading_value)


def internal_any(
    half_span: float = WING_SPAN / 2,
    precision: int = 100,
    **kwargs: np.typing.NDArray | tuple[function, float, float] | tuple[float, float],
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
        if isinstance(val, tuple) and len(val)==2 and isinstance(val[0], float) and isinstance(val[1], float):
            process_type = "point_force"
        elif isinstance(val, np.ndarray):
            distr = ip.interp1d(val[:,0], val[:,1], bounds_error=False, fill_value=0.0)
            process_type = "distr"
        elif callable(val):
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
                    loading[i, 1] += ig.quad(distr, loading[i,0], half_span)[0]
            case _:
                raise ValueError(
                    f"The type of loading '{key}' could not be determined."
                )
    return loading



if __name__ == "__main__":

    """
    ##########################################################################################
    # plt.plot(np.linspace(0,1,30), chord_dist(np.linspace(0,1,30), .2))
    # plt.xlim(0,1)
    # plt.ylim(0,1)
    # plt.show()

    # validation code
    # rd.seed(45)
    # for i in range(10_000):
    #     l = rd.random() * 100
    #     te = rd.random() * l
    #     x = rd.random() * te
    #     weight = 0.5 * FUEL_DENSITY * FUEL_VOLUME_TOTAL * G
    #     weight_ig = ig.quad(lambda x: fuel_weight_distribution(x, l, te), 0, te)
    #     if not abs(weight - weight_ig[0]) < EPSILON:
    #         print(i)
    #         print(weight, weight_ig, weight / weight_ig[0])
    #         break
    # weight_ig = ig.quad(lambda x: wing_weight_distribution(x, l), 0, l)[0]
    # assert abs(WING_MASS * G - weight_ig) < EPSILON
    # print("complete")


    ############################################333##############333
    # # Create the figure and axes
    # fig, ax = plt.subplots(figsize=(10, 4))

    # # Beam dimensions
    # beam_start = 0
    # beam_end = 10

    # # Draw the beam
    # ax.plot(
    #     [beam_start, beam_end],
    #     [0, 0],
    #     color="black",
    #     linewidth=6,
    #     solid_capstyle="butt",
    # )

    # # Draw concentrated forces
    # draw_force(ax, x=3, y_start=2.2, y_end=0.15, label=r"$F_1$")
    # draw_force(ax, x=7, y_start=-2.0, y_end=-0.15, label=r"$F_2$")

    # # Draw applied moments
    # draw_moment(
    #     ax,
    #     x=5,
    #     y=0.8,
    #     radius=0.6,
    #     label=r"$M_1$",
    #     clockwise=False,
    # )

    # draw_moment(
    #     ax,
    #     x=9,
    #     y=0.8,
    #     radius=0.6,
    #     label=r"$M_2$",
    #     clockwise=True,
    # )

    # x_positions = np.linspace(1, 9, 12)

    # for x in x_positions:
    #     draw_force(
    #         ax,
    #         x=x,
    #         y_start=1.5,
    #         y_end=0.15,
    #         # label="",
    #         color="tab:blue",
    #     )

    # # Optional line connecting the arrow tails
    # ax.plot(
    #     [x_positions[0], x_positions[-1]],
    #     [1.5, 1.5],
    #     color="tab:blue",
    #     linewidth=1.5,
    # )

    # ax.text(5, 1.7, r"$w$", ha="center", color="tab:blue")

    # # Plot formatting
    # ax.set_xlim(-1, 11)
    # ax.set_ylim(-2.5, 3)
    # ax.set_aspect("equal")
    # ax.axis("off")

    # fig.tight_layout()
    # plt.show()

    ##############################################################################
    """

    fig, (shear, moment) = plt.subplots(1,2) # 0

    # define loadings
    force_a = (0., -10.)
    force_b = (6., -25.)
    force_2m = (2., 15.)
    moment_3m = (3., -20.)
    force_distr = np.array([[4., 10.], [6., 10.]])

    # calculate internal forces
    shear_loading = internal_any(
        half_span=6.,
        force_a=force_a,
        force_b=force_b,
        force_2m=force_2m,
        force_distr=force_distr,
    )
    moment_loading = internal_any(
        half_span=6.,
        moment_3m=moment_3m,
        moment_shear=np.column_stack((shear_loading[:,0], -shear_loading[:,1]))
    )

    # draw internal forces
    shear.plot(shear_loading[:,0], shear_loading[:,1])
    moment.plot(moment_loading[:,0], moment_loading[:,1])

    # format and show the graphs
    shear.set_xlabel("z")
    moment.set_xlabel("z")
    shear.set_ylabel("internal shear force")
    moment.set_ylabel("interal bending moment")
    plt.show()
