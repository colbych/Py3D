import numpy as np
import matplotlib.pyplot as plt

nx,ny,nz = 256, 256, 64
dx = .05

x = np.arange(nx)*dx + dx/2.
y = np.arange(ny)*dx + dx/2.
z = np.arange(nz)*dx + dx/2.
lx = x[-1] + x[0]
ly = y[-1] + y[0]
lz = z[-1] + z[0]

yy,xx,zz = np.meshgrid(x,y,z)
rr = np.sqrt((xx - lx/2.)**2 + (yy - ly/2.)**2)

psiM = np.exp(-rr**2/8.)

Bx = (roll(psiM,-1,axis=1) - roll(psiM,1,axis=1))/2./dx
By = -(roll(psiM,-1,axis=0) - roll(psiM,1,axis=0))/2./dx
Bm = np.sqrt(Bx**2 + By**2)

levels = linspace(0,1,20)
plt.clf()
plt.contour(x, x, psiM[:,:,0], levels=[.5], colors='k')
#plt.pcolormesh(x, x, psiM)
#plt.pcolormesh(x, x, Bm)
plt.gca().set_aspect('equal')

