
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.interpolate import interp1d
from scipy.integrate import simpson

# XFLR5 export files
folder = Path(__file__).resolve().parent 

# Column indices in the XFLR5 Main Wing table
Y, CHORD, AI, CL, CDI, CM = 0, 1, 2, 3, 5, 7


# Read XFLR5 data for the positive half-wing
def read_xflr5(filename):
    lines = (folder / filename).read_text(encoding="utf-8").splitlines()

    overall_cl = next(
        float(line.split("=", 1)[1])
        for line in lines
        if line.startswith("CL ")
    )

    start = lines.index("Main Wing") + 2
    rows = []

    for line in lines[start:]:
        if not line.strip():
            break
        rows.append([float(value) for value in line.split()])

    data = np.asarray(rows, dtype=float)
    data = data[data[:, Y] > 0]
    data = data[np.argsort(data[:, Y])]

    return data, overall_cl


# Import reference analyses at 0 and 10 degrees
data_0, CL_0 = read_xflr5("MainWing_a=0.00_v=10.00ms.txt")
data_10, CL_10 = read_xflr5("MainWing_a=10.00_v=10.00ms.txt")

# Verify that both files use the same spanwise stations
if data_0.shape != data_10.shape or not np.allclose(
    data_0[:, Y], data_10[:, Y]
):
    raise ValueError("The XFLR5 files have different spanwise stations.")


# Interpolate the local aerodynamic coefficients along the span
span_stations = data_0[:, Y]


def spline(data, column):
    return interp1d(
        span_stations,
        data[:, column],
        kind="cubic",
        bounds_error=False,
        fill_value="extrapolate"
    )


# Chord distribution [m]
c = interp1d(
    span_stations,
    data_0[:, CHORD],
    kind="linear",
    bounds_error=False,
    fill_value="extrapolate"
)

# Local coefficients at 0 degrees
cl_0 = spline(data_0, CL)
cdi_0 = spline(data_0, CDI)
cm_0 = spline(data_0, CM)
ai_0 = spline(data_0, AI)

# Local coefficients at 10 degrees
cl_10 = spline(data_10, CL)
cdi_10 = spline(data_10, CDI)
cm_10 = spline(data_10, CM)
ai_10 = spline(data_10, AI)


# Calculate required whole-wing lift coefficient
def required_cl(load_factor, weight_N, rho, V, wing_area):
    q = 0.5 * rho * V**2
    return load_factor * weight_N / (q * wing_area)


# Interpolation factor between reference analyses
def k_from_cl(CL_target):
    return (CL_target - CL_0) / (CL_10 - CL_0)


# Geometric angle of attack [degrees]
def alpha_deg(CL_target):
    return 10.0 * k_from_cl(CL_target)


# Local lift coefficient [-]
def cl(y_pos, CL_target):
    k = k_from_cl(CL_target)
    return cl_0(y_pos) + k * (cl_10(y_pos) - cl_0(y_pos))


# Local induced angle [degrees]
def ai(y_pos, CL_target):
    k = k_from_cl(CL_target)
    return ai_0(y_pos) + k * (ai_10(y_pos) - ai_0(y_pos))


# Local pitching-moment coefficient about quarter-chord [-]
def cm(y_pos, CL_target):
    k = k_from_cl(CL_target)
    return cm_0(y_pos) + k * (cm_10(y_pos) - cm_0(y_pos))


# Local induced drag coefficient [-]
def cdi(y_pos, CL_target):
    k = k_from_cl(CL_target)

    # Approximate quadratic dependence of induced drag on CL
    # and match the XFLR5 coefficients at 0 and 10 degrees
    delta_cl = cl_10(y_pos) - cl_0(y_pos)
    delta_ai_rad = np.deg2rad(ai_10(y_pos) - ai_0(y_pos))

    return (
        (1 - k) * cdi_0(y_pos)
        + k * cdi_10(y_pos)
        - k * (k - 1) * delta_cl * delta_ai_rad
    )


# Distributed lift [N/m]
def lift(y_pos, rho, V, CL_target):
    q = 0.5 * rho * V**2
    return q * c(y_pos) * cl(y_pos, CL_target)


# Distributed induced drag [N/m]
def drag(y_pos, rho, V, CL_target):
    q = 0.5 * rho * V**2
    return q * c(y_pos) * cdi(y_pos, CL_target)


# Distributed pitching moment [(N m)/m]
def pitching_moment(y_pos, rho, V, CL_target):
    q = 0.5 * rho * V**2
    return q * c(y_pos)**2 * cm(y_pos, CL_target)


# Distributed force normal to the wing chord [N/m]
def normal_force(y_pos, rho, V, CL_target):
    alpha_rad = np.deg2rad(alpha_deg(CL_target))
    return (
        lift(y_pos, rho, V, CL_target) * np.cos(alpha_rad)
        + drag(y_pos, rho, V, CL_target) * np.sin(alpha_rad)
    )


# Flight conditions, checks and plots
if __name__ == "__main__":

    # Flight conditions from WP4.3
    rho = 0.2874             # Air density [kg/m^3]
    V = 200.6                # True airspeed [m/s]
    MTOM = 7990.182          # Maximum take-off mass [kg]
    n = 1.5                  # Load factor [-]
    g = 9.81                 # Gravitational acceleration [m/s^2]
    S = 26.6518587498352     # Total wing area [m^2]

    # Calculate required overall lift coefficient
    W = MTOM * g
    q = 0.5 * rho * V**2
    CL_target = required_cl(n, W, rho, V, S)

    print("Flight conditions")
    print(f"Dynamic pressure = {q:.2f} Pa")
    print(f"Required total lift = {n * W:.2f} N")
    print(f"Required CL = {CL_target:.6f}")
    print(f"Angle of attack = {alpha_deg(CL_target):.3f} deg")

    if not 0 <= k_from_cl(CL_target) <= 1:
        print("WARNING: Extrapolating beyond the 0-10 deg analyses.")

    # Evaluate the functions at one example spanwise position
    y_pos = 4.25
    print(f"\nAerodynamic loads at y = {y_pos:.2f} m")
    print(f"Lift = {lift(y_pos, rho, V, CL_target):.3f} N/m")
    print(f"Drag = {drag(y_pos, rho, V, CL_target):.3f} N/m")
    print(
        f"Pitching moment = "
        f"{pitching_moment(y_pos, rho, V, CL_target):.3f} Nm/m"
    )
    print(
        f"Normal force = "
        f"{normal_force(y_pos, rho, V, CL_target):.3f} N/m"
    )

    # Interpolation plot domain: only the XFLR5 station range
    y_plot = np.linspace(span_stations[0], span_stations[-1], 300)

    L = lift(y_plot, rho, V, CL_target)
    D = drag(y_plot, rho, V, CL_target)
    M = pitching_moment(y_plot, rho, V, CL_target)
    N = normal_force(y_plot, rho, V, CL_target)

    # Plot the four aerodynamic loading distributions
    fig, axes = plt.subplots(2, 2, figsize=(11, 7))

    plots = [
        (L, "Lift distribution", "Lift [N/m]"),
        (D, "Induced drag distribution", "Drag [N/m]"),
        (M, "Pitching moment distribution", "Moment [Nm/m]"),
        (N, "Normal force distribution", "Normal force [N/m]")
    ]

    for ax, (values, plot_title, ylabel) in zip(axes.flat, plots):
        ax.plot(y_plot, values)
        ax.set_title(plot_title)
        ax.set_xlabel("Spanwise position y [m]")
        ax.set_ylabel(ylabel)
        ax.grid(True)

    fig.suptitle(f"Aerodynamic loading")
    plt.tight_layout()

    # Validate overall lift with numerical integration
    # This extrapolates to the physical root and tip of the half-wing
    half_span = 8.55  # [m] Verify against the final wing geometry
    y_check = np.linspace(0.0, half_span, 2000)

    L_integrated = 2 * simpson(
        lift(y_check, rho, V, CL_target), x=y_check
    )

    L_required = n * W
    lift_error_pct = 100 * (L_integrated - L_required) / L_required

    print("\nWhole-wing lift validation")
    print(f"Required lift = {L_required:.2f} N")
    print(f"Integrated lift = {L_integrated:.2f} N")
    print(f"Difference = {lift_error_pct:.2f}%")
    print("Note: check the lift mismatch before using final loads.")

    plt.show()
