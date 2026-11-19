"""
Main file for running strength-limited design simulations
"""

from scipy import integrate as ig
from scipy import interpolate as ip
import numpy as np

import numpy as np

class Line:
    def __init__(self, x, y, L, theta, t):
        self.x = x
        self.y = y
        self.L = L
        self.theta = theta
        self.t = t
    def centroid(self):
        x_cen = self.x
        y_cen = self.y
        
        return x_cen, y_cen
    def Ixx(self):
        # calculation of moment of inertia term
        if theta == (np.pi()/2) or (-np.pi()/2):
            Ixx_s = (1/12)*t*(L**3)
        elif theta == 0:
            Ixx_s = (1/12)*L*(t**3)
        else:
            Ixx_s = (1/12)*t*(L**3)*((np.sin(theta))**2)
        # calculation of steiner terms
        Area = L*t
        Ixx_st = Area * y ** 2

        return Ixx
    def Iyy(self):

class Point:
    def __init__(self, x, y, A):
        self.x = x
        self.y = y
    def centroid(self):
    def Ixx(self):
    def Iyy(self):