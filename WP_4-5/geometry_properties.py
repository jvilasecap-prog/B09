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
        if abs(self.theta) == np.pi/2:
            Ixx_s = (1/12)*self.t*(self.L**3)
        elif self.theta == 0:
            Ixx_s = (1/12)*self.L*(self.t**3)
        else:
            Ixx_s = (1/12)*self.t*(self.L**3)*((np.sin(self.theta))**2)
        # calculation of steiner terms
        Area = self.L*self.t
        Ixx_st = Area * abs(self.y-D_CEN[1]) ** 2
        # summing both components
        Ixx = Ixx_s + Ixx_s
        return Ixx
    def Iyy(self, D_CEN):
        # calculation of moment of inertia term
        if abs(self.theta) == np.pi/2:
            Iyy_s = (1/12)*self.t*(self.L**3)
        elif self.theta == 0:
            Iyy_s = (1/12)*self.L*(self.t**3)
        else:
            Iyy_s = (1/12)*self.t*(self.L**3)*((np.cos(self.theta))**2)
        # calculation of steiner terms
        Area = self.L*self.t
        Iyy_st = Area * abs(self.x-D_CEN[0]) ** 2
        # summing both components
        Iyy = Iyy_s + Iyy_s
        return Iyy

class Point:
    def __init__(self, x, y, A):
        self.x = x
        self.y = y
        self.A = A
    def Ixx(self, D_CEN):
        I_xx = self.area * (abs(self.y-D_CEN[1])) ** 2
        return I_xx
    def Iyy(self, D_CEN):
        I_yy = self.area * (abs(self.x-D_CEN[0])) ** 2
        return I_yy