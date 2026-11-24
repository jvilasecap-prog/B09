import numpy as np
from geometry_properties import Line, Point

#Get the centroid, Ixx and Iyy of a geometry
#The convention is x left y down, the wing box is a trapezoid with the following points
#A - B
#D - C

#Chord is a function of the c

#MAKE A CODE TO GET THE LENGTH, X, Y, THETA POSITIONS OF AB, BC, CD, AD
#SO THAT IT CAN BE USED BY THE CLASSES

#AB, BC, CD, AD are arrays of 2 numbers (x, y)
def get_L(A, B):
  L = np.sqrt((B[0]-A[0])**2+(B[1]-A[1])**2)
  return L
def get_theta(A, B):
  if (B[0]- A[0]) == 0:
     if (B[1] - A[1]) > 0:
        theta = np.pi/2
     if (B[1] - A[1]) < 0:
        theta = -np.pi/2
  else:
    theta = np.atan((B[1]-A[1])/(B[0]-A[0]))
  return theta

def get1_L_value(A, B, t):
    L = get_L(A, B)
    theta = get_theta(A, B)
    x_mid = (A[0] + B[0]) / 2
    y_mid = (A[1] + B[1]) / 2
    line_added = Line(x_mid, y_mid, L, theta, t)
    return line_added

def get_L_values(A, B, C, D, t):
  Line_AB = get1_L_value(A, B, t)
  Line_BC = get1_L_value(B, C, t)
  Line_CD = get1_L_value(C, D, t)
  Line_DA = get1_L_value(D, A, t)
  return Line_AB, Line_BC, Line_CD, Line_DA

# there will be an array size [N, 3] where N is the number of area points and [3] is [x, y, A]
def get_A_values(A_array):
    return [Point(float(x), float(y), float(a)) for x, y, a in A_array]
  
def global_Centroid(Lines, Points):
    line_area = sum(line.L * line.t for line in Lines)
    pt_area = sum(pt.A for pt in Points)
    total_area = line_area + pt_area

    Qx = sum(line.x * line.L * line.t for line in Lines) + sum(pt.x * pt.A for pt in Points)
    Qy = sum(line.y * line.L * line.t for line in Lines) + sum(pt.y * pt.A for pt in Points)

    return np.array([Qx / total_area, Qy / total_area])
 
def global_I(Lines, Points, x_cen):
  Ixx = sum(Line.Ixx(x_cen) for Line in Lines) + sum(Point.Ixx(x_cen) for Point in Points)
  Iyy = sum(Line.Iyy(x_cen) for Line in Lines) + sum(Point.Iyy(x_cen) for Point in Points)
  return Ixx, Iyy


  
  