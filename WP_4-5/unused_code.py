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

