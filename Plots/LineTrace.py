
# Fourth test of field line tracing in 2D

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

# subsampling is necessary if plotting over Jz
def SubSample(Par_Orig, rate):
    print('subsampling....')
    Par_Sub = np.zeros((8192/rate, 8192/rate))
    for i in range(0,8192/rate):
        for n in range(0,8192/rate):
            Par_Sub[i,n] = Par_Orig[i*rate,n*rate]
    return Par_Sub

# SS_Rate may be set to 1 to avoid subsampling without disturbing code
SS_Rate = 8
# subsampling calls if necessary
d['jz'] = SubSample(d['jz'], SS_Rate)
d['bx'] = SubSample(d['bx'], SS_Rate)
d['by'] = SubSample(d['by'], SS_Rate)
#d['bz'] = SubSample(d['bz'], SS_Rate)

# the magnetic field components
By = d['bx']
Bx = d['by']

# line tracing method
def Line(InitX, InitY): 
    
    # maximum number of steps along line if it does not hit an edge
    MaxSteps = 1000000
    
    # slope variable to hold Bx/|B| at point (X,Y)
    DeltaX = 0
    # and By/|B|
    DeltaY = 0
    # variables to hold the slopes from point (X-1,Y-1) at the previous step
    DeltaX_Store = 0
    DeltaY_Store = 0
    
    # arrays to hold X and Y coordinates of field line points
    Line_X = np.zeros(MaxSteps)
    Line_Y = np.zeros(MaxSteps)
    
    # initial point for a field line
    X = InitX
    Y = InitY
    
    # real number of steps along line (incase it hits an edge)
    Steps = 0
    
    # position holders to avoid recalculating nearest neighbors if still
    # within the same grid space between steps
    iStore = -1
    jStore = -1
    
    # loop to step forward line from initial point
    for step in range(0,MaxSteps):
        # trimming line for handling hitting edge
        if X <= 1 or Y <= 1 or Y >= 8191/SS_Rate or X >= 8191/SS_Rate:
            Line_X = Line_X[:Steps]
            Line_Y = Line_Y[:Steps]
            break
        # counter for trimming if edge is hit
        Steps = Steps + 1
        # just a status update, sometimes it takes a long time
        if (Steps % 1000000) == 0:
            print('Step = ',Steps)
        
        # adding points (X,Y) to field line
        Line_X[step] = X
        Line_Y[step] = Y
        
        # differential step towards next point on field line
        # should this change with subsampling?
        dx = .0025
        
        # distance between field line point (X,Y) and lower left data grid 
        #point (i,j), such that Wx = X - i, Wy = Y - j
        Wx = X % 1
        Wy = Y % 1
        
        # corresponding closest lower left data grid point (i,j) to field line
        # point (X,Y)
        i = int(X - Wx)
        j = int(Y - Wy)
        
        # checking to see if new grid space has been entered
        if i != iStore or j != jStore:
            # identifying Bx, By at closest 4 data grid points
            Bx_ij = Bx[i,j]
            Bx_i1j = Bx[i+1,j]
            Bx_ij1 = Bx[i,j+1]
            Bx_i1j1 = Bx[i+1,j+1]
        
            By_ij = By[i,j]
            By_i1j = By[i+1,j]
            By_ij1 = By[i,j+1]
            By_i1j1 = By[i+1,j+1]
            
            # stores coordinates of lower left neighbor to test if neighbors
            # have changed between steps so they don't need to be recalculated
            # saves considerable run time
            iStore = i
            jStore = j
        
        # finding average Bx, By and Bm from 4 nearest neighboring data points
        # (i,j), (i+1,j), (i,j+1) and (i+1,j+1) at field line point (X,Y)
        B_Wx = (1-Wx)*(1-Wy)*Bx_ij + (1-Wx)*Wy*Bx_ij1 + Wx*(1-Wy)*Bx_i1j + Wx*Wy*Bx_i1j1
        B_Wy = (1-Wx)*(1-Wy)*By_ij + (1-Wx)*Wy*By_ij1 + Wx*(1-Wy)*By_i1j + Wx*Wy*By_i1j1
        B_Wm = np.sqrt(B_Wx**2 + B_Wy**2)
        
        # finding next point on field line
        # first calculate slope of B at (X,Y)
        DeltaX = ((dx*B_Wx/B_Wm))
        DeltaY = ((dx*B_Wy/B_Wm))
        if Steps > 1:
            # then find next point from slope at (X,Y) -> (Delta) and 
            # (X-1,Y-1) -> (Delta_Store)
            Y = Y + DeltaY + .5*(DeltaY - DeltaY_Store)
            X = X + DeltaX + .5*(DeltaX - DeltaX_Store)
        # on first step we have no slope from (X-1,Y-1)
        else:
            Y = Y + DeltaY
            X = X + DeltaX
        # storing slope of (X,Y) for next iteration
        DeltaY_Store = DeltaY
        DeltaX_Store = DeltaX
        
    # returns line as two arrays of X and Y points
    return Line_X, Line_Y

# Plotting single field line
# Contour = Line(4000,4000)
#ax.plot(Contour[0],Contour[1])

# dictionary is used to conveiniently calculate 64 lines starting at points in
# an evenly spaced grid across the 2-D field
Lines = dict()
for t in range(1,9):
    for q in range(1,9):
        print('Tracing Line # ',((t-1)*8)+q,' of 64')
        Lines[((t-1)*8)+q] = Line(((1000/SS_Rate)*t)-(404/SS_Rate),((1000/SS_Rate)*q)-(404/SS_Rate))


# plotting the field lines
for i in range(1,65):
    # plotting Jz if desired
    ax.pcolormesh(d['jz'].T)
    ax.plot(Lines[i][0],Lines[i][1])

plt.show()
