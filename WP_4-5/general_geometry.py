import numpy as np
from geometry_properties import Line, Point

#Get the centroid, Ixx and Iyy of a geometry
#The convention is x left y down, the wing box is a trapezoid with the following points
#A - B
#D - C
LineAB = Line(0.5, 0, 1, 0, 0.3)