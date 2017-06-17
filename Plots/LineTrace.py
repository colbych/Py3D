
# Second test of field line tracing in 2D

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


import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
from Py3D.sub import load_movie

# the figure
fig1 = plt.figure(1)
fig1.set_size_inches(10,10, forward = True)
ax = fig1.add_subplot(111)  #, projection = '3d')

# loading data
d = load_movie()

def SubSample(Par_Orig, rate):
    Par_Sub = np.zeros((8192/rate, 8192/rate))
    for i in range(0,8192/rate):
        for n in range(0,8192/rate):
            Par_Sub[i,n] = Par_Orig[i*rate,n*rate]
    return Par_Sub

SS_Rate = 1
d['jiz'] = SubSample(d['jiz'], SS_Rate)
d['jez'] = SubSample(d['jez'], SS_Rate)
d['bx'] = SubSample(d['bx'], SS_Rate)
d['by'] = SubSample(d['by'], SS_Rate)
d['bz'] = SubSample(d['bz'], SS_Rate)

Bx = d['bx']
By = d['by']
Bm = np.sqrt(Bx**2 + By**2)

YStore = 0
def Line(InitX, InitY):  
    # number of steps along line
    MaxSteps = 8000
    
    # the actual field line
    Line_X = np.zeros(MaxSteps)
    Line_Y = np.zeros(MaxSteps)
    
    # initial point for a field line
    X = InitX
    Y = InitY
    # real number of steps along line (incase it hits an edge)
    Steps = 0
    # loop to step forward line from initial point
    for i in range(0,MaxSteps):
        # trimming line for handling hitting edge
        if X <= 0 or Y <= 0 or Y >= 8192/SS_Rate or X >= 8192/SS_Rate:
            Line_X = Line_X[:Steps]
            Line_Y = Line_Y[:Steps]
            break
        Steps = Steps + 1
        # adding points to line
        Line_X[i] = X
        Line_Y[i] = Y
        # points to skip towards next point on line
        P = 10
        # finding next point on line
        YStore = Y
        Y = Y + ((P*By[X,Y]/Bm[X,Y]) - (P*By[X,Y]/Bm[X,Y]) % 1)
        X = X + ((P*Bx[X,YStore]/Bm[X,YStore]) - (P*Bx[X,YStore]/Bm[X,YStore]) % 1)
        
    return Line_X, Line_Y

Lines = dict()
for t in range(1,9):
    for q in range(1,9):
        Lines[((t-1)*8)+q] = Line(((1000/SS_Rate)*t)-(404/SS_Rate),((1000/SS_Rate)*q)-(404/SS_Rate))

for i in range(1,65):
    #ax.pcolormesh(d['jiz']+d['jez'])
    ax.plot(Lines[i][0],Lines[i][1])

plt.show()
