
# first test of field line tracing in 2D

# www2.warwick.ac.uk
#   - good resource explaining definition of magnetic field lines
#   - field lines satisfy coupled differential equations
#     dx/ds = Bx/B , dy/ds = By/B , dz/ds = Bz/B
#     where B is magnitude of magnetic field and x, y, and z are components,
#     and ds is a differential step along the field line
#   - these equations describe the slope / direction of the magnetic field 
#     which in unit vector form is < Bx, By, Bz > / B, corresponding to the
#     slope of the field line < dx, dy, dz > / ds
#   - Initial spacing of field lines will be decided by choosing set of
#     evenly spaced initial points
#   - these points can then be stepped forward by checking the slope at each point

# http://www.math.poly.edu/courses/ma2132/Notes/MA2132NumericalMethods.pdf
#   - solid resource for Euler's midpoint method
# https://www.saylor.org/site/wp-content/uploads/2011/06/MA221-6.1.pdf
#   - solid resource for Verlet's method

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
from Py3D.sub import load_movie

# the figure
fig1 = plt.figure(1)
ax = fig1.add_subplot(111)  #, projection = '3d')

# loading data
d = load_movie()

SS_Rate = 8

def SubSample(Par_Orig, rate):
    Par_Sub = np.zeros((8192/rate, 8192/rate))
    for i in range(0,8192/rate):
        for n in range(0,8192/rate):
            Par_Sub[i,n] = Par_Orig[i*rate,n*rate]
    return Par_Sub


d['jiz'] = SubSample(d['jiz'], SS_Rate)
d['bx'] = SubSample(d['bx'], SS_Rate)
d['by'] = SubSample(d['by'], SS_Rate)
d['bz'] = SubSample(d['bz'], SS_Rate)

# creating magnetic field component arrays
Comp_X = np.zeros((8192/SS_Rate,8192/SS_Rate))
Comp_Y = np.zeros((8192/SS_Rate,8192/SS_Rate))
#Comp_Z = np.zeros((8192/SS_Rate,8192/SS_Rate))

# Magnitude of magnetic field
Mb  = np.sqrt(d['bx']**2 + d['by']**2)  #+ d['bz']**2)

# normalizing the components
Comp_X = d['bx']/Mb
Comp_Y = d['by']/Mb
#Comp_Z = d['bz']/Mb

def Line(InitX, InitY):  
    # number of steps along line
    # need to figure out way to handle reaching the egde of data
    LineSteps = 10000
    
    # the actual field line
    Line_X = np.zeros(LineSteps)
    Line_Y = np.zeros(LineSteps)
    
    # initial point for a field line
    X = InitX
    Y = InitY
    # initial data point of field line
    Line_X[0] = X
    Line_Y[0] = Y
    # loop to step forward line from initial point
    for i in range(1,LineSteps):
        if X > (8192/SS_Rate)-1:
            X = (8192/SS_Rate)-1
        if Y > (8192/SS_Rate)-1:
            Y = (8192/SS_Rate)-1
        if X < 0:
            X = 0
        if Y < 0:
            Y = 0
        # stepping to next point based on magnetic field slope
        # angle of magnetic field at point (X,Y)
        Angle = np.rad2deg(np.arctan2(Comp_Y[X, Y] , Comp_X[X, Y]))
        # to handle output format of np.acrtan2(x,y)
        if Angle < 0:
            Angle = Angle + 360
        # first check quadrant for easy computation of slope / angle
        if 0 <= Angle <= 90:
            # if between 30 and 60 degrees
            if 30 <= Angle <= 60:
                # go to upper right corner
                X = X + 1
                Y = Y + 1
            elif Angle < 30:
                # go right
                X = X +1
            elif Angle > 60:
                # go up
                Y = Y + 1
        elif 90 <= Angle <= 180:
            if 120 <= Angle <= 150:
                # go to upper left corner
                X = X - 1
                Y = Y + 1
            elif Angle < 120:
                # go up
                Y = Y + 1
            elif Angle > 150:
                # go left
                X = X - 1
        elif 180 <= Angle <= 270:
            if 210 <= Angle <= 240:
                # go to lower left corner
                X = X - 1
                Y = Y - 1
            elif Angle < 210:
                # go left
                X = X - 1
            elif Angle > 240:
                # go down
                Y = Y - 1
        elif 270 <= Angle <= 360:
            if 300 <= Angle <= 330:
                # go to lower right corner
                X = X + 1
                Y = Y - 1
            elif Angle < 300:
                # go down
                Y = Y - 1
            elif Angle > 330:
                # go right
                X = X + 1
        Line_X[i] = X
        Line_Y[i] = Y
    return Line_X, Line_Y

Lines = dict()
for t in range(1,9):
    for q in range(1,9):
        Lines[((t-1)*8)+q] = Line(((1000/SS_Rate)*t)-(404/SS_Rate),((1000/SS_Rate)*q)-(404/SS_Rate))

for i in range(1,65):
    ax.pcolormesh(d['jiz'])
    ax.plot(Lines[i][0],Lines[i][1])

plt.show()
