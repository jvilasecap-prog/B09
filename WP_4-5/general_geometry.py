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
  theta = np.atan((B[1]-A[1])/(B[0]-A[0]))

def get_L_values(A, B, C, D, t):
  L_AB = get_L(A, B)
  L_BC = get_L(B, C)
  L_CD = get_L(C, D)
  L_DA = get_L(D, A)
  theta_AB = get_theta(A, B)
  theta_BC = get_theta(B, C)
  theta_CD = get_theta(C, D)
  theta_DA = get_theta(D, A)
  Line_AB = Line((A[0] + B[0])/2, A[1] + B[1]/2, L_AB, theta_AB, t)
  Line_BC = Line((B[0] + C[0])/2, B[1] + C[1]/2, L_BC, theta_AB, t)
  Line_CD = Line((C[0] + D[0])/2, C[1] + D[1]/2, L_CD, theta_AB, t)
  Line_DA = Line((D[0] + A[0])/2, D[1] + B[1]/2, L_DA, theta_AB, t)

  return Line_AB, Line_BC, Line_CD, Line_DA
def global_Centroid(Lines):
  n = len(Lines)
  cx = sum(Line.x for Line in Lines) / n
  cy = sum(Line.y for Line in Lines) / n

  return cx,cy
# there will be an array size [N, 3] where N is the number of area points and [3] is [x, y, A] 
def global_Ixx(Lines):
  n = len(Lines)
  Ixx = sum(Ixx)
def get_A_values(A_array):
  
  