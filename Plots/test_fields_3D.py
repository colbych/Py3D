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

Bx = (np.roll(psiM,-1,axis=1) - np.roll(psiM,1,axis=1))/2./dx
By = -(np.roll(psiM,-1,axis=0) - np.roll(psiM,1,axis=0))/2./dx
#Bz = psiM
Bz = 0.*zz + 1.

Bm = np.sqrt(Bx**2 + By**2 + Bz**2)

def dif(i,f):
    return (np.roll(f,-1,axis=i) - np.roll(f,1,axis=i))/2./dx

# You can comment the divergence and curl parts out if you want.

# Calculate the Divergence, make sure it is 0
DB = dif(0, Bx) + dif(1, By) + dif(2, Bz)

# Calculate the Curl, this gives J
CBx = dif(1, Bz) - dif(2, By) 
CBy = dif(2, Bx) - dif(0, Bz) 
CBz = dif(0, By) - dif(1, Bx) 

plt.clf()
# This is no longer the correct field line in the strict sence
plt.contour(x, y, psiM[:,:,0], levels=[.125,.25,.5,1.,2.], colors='k')
plt.pcolormesh(x, y, Bm[:,:,0])
plt.gca().set_aspect('equal')
plt.show()
