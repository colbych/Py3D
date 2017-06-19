
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
SS_Rate = 1
# subsampling calls if necessary
#d['jz'] = SubSample(d['jz'], SS_Rate)
#d['bx'] = SubSample(d['bx'], SS_Rate)
#d['by'] = SubSample(d['by'], SS_Rate)
#d['bz'] = SubSample(d['bz'], SS_Rate)

# the magnetic field components
Bx = d['bx']
By = d['by']

# line tracing method
def Line(InitX, InitY): 
    
    # maximum number of steps along line if it does not hit an edge
    MaxSteps = 50000000
    
    # initial slope parameters
    DeltaX = 0
    DeltaY = 0
    DeltaX_Store = 0
    DeltaY_Store = 0
    
    # the actual field line
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
        if X <= 0 or Y <= 0 or Y >= 8192/SS_Rate or X >= 8192/SS_Rate:
            Line_X = Line_X[:Steps]
            Line_Y = Line_Y[:Steps]
            break
        # counter for hitting an edge
        Steps = Steps + 1
        if (Steps % 1000000) == 0:
            print('Step = ',Steps)
        
        # adding points (X,Y) to field line
        Line_X[step] = X
        Line_Y[step] = Y
        
        # differential step towards next point on field line
        # should this change with subsampling?
        d = .00125
        
        # distance between field line point (X,Y) and data grid point (i,j)
        Wx = X % 1
        Wy = Y % 1
        
        # corresponding closest lower left data grid point (i,j)
        i = X - Wx
        j = Y - Wy
        
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
            
            iStore = i
            jStore = j
        
        # finding average Bx, By and Bm at (X,Y)
        B_Wx = (1-Wx)*(1-Wy)*Bx_ij + (1-Wx)*Wy*Bx_ij1 + Wx*(1-Wy)*Bx_i1j + Wx*Wy*Bx_i1j1
        B_Wy = (1-Wx)*(1-Wy)*By_ij + (1-Wx)*Wy*By_ij1 + Wx*(1-Wy)*By_i1j + Wx*Wy*By_i1j1
        B_Wm = np.sqrt(B_Wx**2 + B_Wy**2)
        
        # finding next point on line
        DeltaX = ((d*B_Wx/B_Wm))
        DeltaY = ((d*B_Wy/B_Wm))
        if Steps > 1:
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

Contour = Line(4000,4000)
# dictionary is used to conveiniently calculate 64 lines starting at points in
# an evenly spaced grid across the 2-D field
#Lines = dict()
#for t in range(1,9):
#    for q in range(1,9):
#        print('Tracing Line # ',((t-1)*8)+q,' of 64')
#        Lines[((t-1)*8)+q] = Line(((1000/SS_Rate)*t)-(404/SS_Rate),((1000/SS_Rate)*q)-(404/SS_Rate))


# plotting the field lines
#for i in range(1,65):
    # plotting Jz if desired
ax.pcolormesh(d['jz'])
    #ax.plot(Lines[i][0],Lines[i][1])
ax.plot(Contour[0],Contour[1])

plt.show()
