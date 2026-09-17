import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import flaps

wing = {
    "Chord": 1.5,             # Mean aerodynamic chord [m]
    "CL_alpha_clean": 0.1,    # Clean lift curve slope [1/deg]
    "CL_max": 1.6,           # Clean maximum lift coefficient [-]
    "a0L": -2.0,              # Zero-lift angle of attack [deg]
    "choices": {
        "flap_type": "single slotted fowler flap"
    }
}

flaps.flaps_calculation(wing)