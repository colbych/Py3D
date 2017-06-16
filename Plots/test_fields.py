import numpy as np
import matplotlib.pyplot as plt
dx = .05
x = np.arange(256)*dx + dx/2.
lx = x[-1] + x[0]

yy,xx = np.meshgrid(x,x)
rr = np.sqrt((xx - lx/2.)**2 + (yy - lx/2.)**2)

psi = np.exp(-rr**2/8.)

i0,j0 = 128, 88
x0,y0 = x[i0], y[j0]
psi0 = psi[i0, j0]

Bx = (roll(psi,-1,axis=1) - roll(psi,1,axis=1))/2./dx
By = -(roll(psi,-1,axis=0) - roll(psi,1,axis=0))/2./dx
Bm = np.sqrt(Bx**2 + By**2)

levels = linspace(0,1,20)
plt.clf()
subplot(211)
plt.contour(x, x, psi, levels=[psi0], colors='k')
plt.gca().set_aspect('equal')

subplot(212)
#plt.pcolormesh(x, x, psi)
plt.pcolormesh(x, x, Bm)
plt.contour(x, x, psi, levels=levels, colors='w', linewidths=.5)
plt.contour(x, x, psi, levels=levels, colors='k', linewidths=.5, linestyles='dashed')
plt.gca().set_aspect('equal')

