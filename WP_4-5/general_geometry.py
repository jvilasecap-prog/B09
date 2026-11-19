import numpy as np
from geometry_properties import Line, Point

#Get the centroid, Ixx and Iyy of a geometry
#The convention is x left y down, the wing box is a trapezoid with the following points
#A - B
#D - C

#Chord is a function of the c

#MAKE A CODE TO GET THE LENGTH, X, Y, THETA POSITIONS OF AB, BC, CD, AD
#SO THAT IT CAN BE USED BY THE CLASSES
def _get_chord_values(c, b1, b2, ):

LineAB = Line(0.5, 0, 1, 0, 0.3)
