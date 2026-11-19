"""
Main file for running strength-limited design simulations
"""

from scipy import integrate as ig
from scipy import interpolate as ip
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
    def Ixx(self, D_CEN):
        # calculation of moment of inertia term
        if theta == (np.pi()/2) or (-np.pi()/2):
            Ixx_s = (1/12)*t*(L**3)
        elif theta == 0:
            Ixx_s = (1/12)*L*(t**3)
        else:
            Ixx_s = (1/12)*t*(L**3)*((np.sin(theta))**2)
        # calculation of steiner terms
        Area = L*t
        Ixx_st = Area * abs(y-D_CEN[1]) ** 2
        # summing both components
        Ixx = Ixx_s + Ixx_s
        return Ixx
    def Iyy(self, D_CEN):
        # calculation of moment of inertia term
        if theta == (np.pi()/2) or (-np.pi()/2):
            Iyy_s = (1/12)*t*(L**3)
        elif theta == 0:
            Iyy_s = (1/12)*L*(t**3)
        else:
            Iyy_s = (1/12)*t*(L**3)*((np.cos(theta))**2)
        # calculation of steiner terms
        Area = L*t
        Iyy_st = Area * abs(x-D_CEN[0]) ** 2
        # summing both components
        Iyy = Iyy_s + Iyy_s
        return Iyy

class Point:
    def __init__(self, x, y, A):
        self.x = x
        self.y = y
    def centroid(self):
    def Ixx(self):
    def Iyy(self):