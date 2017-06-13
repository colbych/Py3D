
# this module codes the test particle run base

import os
import numpy as np
import scipy.ndimage as ndimage

# now for linking C to python
from numpy.ctypeslib import ndpointer
import ctypes

# for plotting routines
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class TPRun:
    """
    Test Particle run
    """

    # This is called here so realpath is where the .so file is
    _pathlibdir   = os.path.dirname(os.path.realpath(__file__))

    #==========================================================
    #==========================================================
    def __init__(self,
                 CR,
                 npart=1,
                 charge=-1.,
                 mass=0.04,
                 run=000,
                 tstart=0,
                 tend=5.,
                 t0=.0,
                 r0=[130.,30.0],
                 dr0=[0.1,0.1],
                 v0=[0.,0.,0.],
                 dv0=[0.1,0.1],
                 dt=.001,
                 loading     = 'randu',     #randu,randn,copy
                 fieldinterp = False):

        """ constructor of the testparticle run object

        @param CR          : The IDL restore file with all the data
        @param npart       : # of particles
        @param charge      : charge of the particles
        @param mass        : mass of the particles
        @param run         : run from which the fiels will be read
        @param tstart      : starting time of the integration
        @param tend        : end time of the integration
        @param t0          : time at which positions/velocities are given
        @param r0          : initial position
        @param dr0         : initial spatial deviation
        @param v0          : initial velocity
        @param dv0         : initial velocity deviation
        @param dt          : time step
        @param loading     : either 'randu', 'randu_mag', 'randn', or 'copy'
        @param fieldinterp : delfault(False) inteprolate fields in time or not

        if loading == 'randu' particles will be loaded randomly in a rectangle
        of size dr0
        if 'randn' is chosen, they will be loaded in a gaussian of spatial std=dr0
        if 'copy' is chosen, r0,v0 have to eb arrays of size (3,npart), dr0 and dv0
        are then disregarded


        @return: a TPRun object

        Creation : 2013-05-01 11:17:42.987369

        """
        self._CR        = CR
        self._npart     = npart
        self._charge    = charge
        self._mass      = mass
        self._r0        = None
        self._v0        = None
        self._r0        = np.array([[r0[0]],[r0[1]],[0.]])
        self._v0        = np.array([[v0[0]],[v0[1]],[v0[2]]])
        self._is3D      = (np.ndim(self._CR['bzav']) == 3)

        # keep the user values for plotting etc.

        #TODO attention si loading==user
        # r0,v0 sont des tableaux et dr0,dv0
        # ne doivent pas etre utilises !
        #- 2013-05-20 08:06:57.972466
        self._r0u       = r0
        self._v0u       = v0
        self._dr0u      = dr0
        self._dv0u      = dv0

        self.dt         = dt
        self.tstart     = tstart
        self.tend       = tend

        self.t0         = t0
        self._nt        = int((tend-tstart)/dt) + 1


        # checks the time interval is well defined
        #----------------------------------------------------------------------
        if t0 < self.tstart or t0 > self.tend:
            print 'time (%5.3f) should be between tstart(%5.3f)\
                   and tend(%5.3f) (both included)' \
                   % (t0 ,self.tstart,self.tend)

            return None
        #----------------------------------------------------------------------
        # that is the selection time index.
        self._it0 = int((self.t0   - self.tstart)/self.dt) + 1


        # loading method
        self._loading = loading
        if loading.lower() == 'randu':
            self.load_randu(r0,dr0,v0,dv0)

        elif loading.lower() == 'randu_mag':
            self.load_randu_mag(r0,dr0,v0,dv0)

        elif loading.lower() == 'randn':
            self.load_randn(r0,dr0,v0,dv0)

        #c# elif loading.lower() == 'user':
        #c#     self._r0 = np.zeros((3,self._npart))
        #c#     self._r0[0,:] = r0[0,:]
        #c#     self._r0[1,:] = r0[1,:]
        #c#     self._v0 = np.zeros((3,self._npart))
        #c#     self._v0 = v0

        #c# else:
        #c#     print 'Ptest : warning, no loading method specified'


        # position and velocity arrays

        self.r   = np.zeros((3, self._npart, self._nt),
                                    order='FORTRAN',
                                    dtype=np.float64)

        self.v   = np.zeros((3, self._npart, self._nt),
                                    order='FORTRAN',
                                    dtype=np.float64)


        # electric and magnetic field seen by each particle
        self.Ep   = np.zeros((3, self._npart, self._nt),
                             order = 'FORTRAN',
                             dtype=np.float64)

        self.Bp   = np.zeros((3, self._npart, self._nt),
                             order = 'FORTRAN',
                             dtype = np.float64)


        # do we interpolate fields in time ?
        self._fldinterp = fieldinterp

        # the PIC run from which he get the fields
        self._run       = run


        # if we interpolate fields in time
        # we actually need to read all the files between the time
        # interval
#cc# Right now we are just assuming that we have an idl file correctly loaded in
#cc# This would be a good point to add in stuf that would use the p3d run object
        if self._fldinterp == False:
            #self._E = np.concatenate(([CR['exav']],[CR['eyav']],[CR['ezav']]), axis=0)
            #self._B = np.concatenate(([CR['bxav']],[CR['byav']],[CR['bzav']]), axis=0)

# This is iffy because we have some fields structured as (y, x) 
# and some fields with (x, y, z)
            if self._is3D:
                self._E = np.zeros((3,self._CR['bzav'].shape[0],
                                      self._CR['bzav'].shape[1],
                                      self._CR['bzav'].shape[2],),
                                      "float32",order='FORTRAN')

                self._B = np.zeros((3,self._CR['bzav'].shape[0],
                                      self._CR['bzav'].shape[1],
                                      self._CR['bzav'].shape[2],),
                                      "float32",order='FORTRAN')
                self._E[0] = self._CR['exav']
                self._E[1] = self._CR['eyav']
                self._E[2] = self._CR['ezav']
                self._B[0] = self._CR['bxav']
                self._B[1] = self._CR['byav']
                self._B[2] = self._CR['bzav']

            else:
                self._E = np.zeros((3,self._CR['bzav'].shape[1],
                                      self._CR['bzav'].shape[0]),
                                      "float32",order='FORTRAN')
                self._B = np.zeros((3,self._CR['bzav'].shape[1],
                                      self._CR['bzav'].shape[0]),
                                      "float32",order='FORTRAN')

                self._E[0,:,:] = self._CR['exav'].transpose()
                self._E[1,:,:] = self._CR['eyav'].transpose()
                self._E[2,:,:] = self._CR['ezav'].transpose()
                self._B[0,:,:] = self._CR['bxav'].transpose()
                self._B[1,:,:] = self._CR['byav'].transpose()
                self._B[2,:,:] = self._CR['bzav'].transpose()
            #c# self._E = np.array([transpose(CR['exav']),transpose(CR['eyav']),transpose(CR['ezav'])])
            #c# self._B = np.array([transpose(CR['bxav']),transpose(CR['byav']),transpose(CR['bzav'])])

            # smooth fields may help... B is fine, E is noisy !!
            sig = 0.
            for c in range(3):
                self._E[c,:,:] = ndimage.gaussian_filter(self._E[c,:,:],
                                                         sigma=sig, #used to be 6
                                                         order=0)

            # checks which component is the out of plane
            #c# if self._run.outofplane == 1:
            #c#     Ey = self._E[1,:,:].copy()
            #c#     Ez = self._E[2,:,:].copy()
            #c#     self._E[1,:,:] = np.copy(Ez)
            #c#     self._E[2,:,:] = np.copy(-Ey)

            #c#     By = self._B[1,:,:].copy()
            #c#     Bz = self._B[2,:,:].copy()
            #c#     self._B[1,:,:] = np.copy(Bz)
            #c#     self._B[2,:,:] = np.copy(-By)



        # in case we want to debug, use analytic fields
        debug = False

        if debug == True:
            self._B[0,:,:]  = 0.#1. + np.random.randn(self._B.shape[1],self._B.shape[2])*1e-3
            self._B[1,:,:]  = 0.#np.random.randn(self._B.shape[1],self._B.shape[2])*1e-3
            self._B[2,:,:]  = 1.#np.random.randn(self._B.shape[1],self._B.shape[2])*1e-3
            self._E[0,:,:]  = 0.#np.random.randn(self._B.shape[1],self._B.shape[2])*1e-2
            self._E[1,:,:]  = 1.#np.random.randn(self._B.shape[1],self._B.shape[2])*1e-2
            self._E[2,:,:]  = 0.#np.random.randn(self._B.shape[1],self._B.shape[2])*1e-2


        # linking to the C library
        pathlib   = os.path.join(self._pathlibdir,'functions.so')
        self._lib = ctypes.cdll.LoadLibrary(pathlib)
    #==========================================================

    #==========================================================
    def move(self):
        """moves all the particles

        @return: @todo

        Exemple  :

        Creation : 2013-05-01 14:28:08.724842

        """

        #c# xr = self._run.GetCoord(axis=0)
        #c# yr = self._run.GetCoord(axis=1)

        #c# xmin = xr[0]
        #c# ymin = yr[0]

        # setup the initial condition

        self.r[:,:,self._it0] = self._r0
        self.v[:,:,self._it0] = self._v0

        # call super-fast cython functions now haha

        # the following function will move all the particles
        # for all time steps
        if self._is3D:
            print 'calling _moveall3D'
            self._moveall3D()
            print 'calling _pfields3D'
            self._pfields3D()

        else:
            print 'calling _moveall'
            self._moveall()
            print 'calling _pfields'
            self._pfields()

    #==========================================================




    #==========================================================
    #==========================================================
    def _moveall(self):
        """move all particles backward and forward as necessary

        Exemple  :

        Creation : 2013-05-04 15:40:32.257893

        """
        func          = self._lib.moveall
        func.restype  = None
        func.argtypes = [ndpointer(ctypes.c_double),    # pos
                         ndpointer(ctypes.c_double),    # vel
                         ndpointer(ctypes.c_float),     # E
                         ndpointer(ctypes.c_float),     # B
                         ctypes.c_uint,                 # nx
                         ctypes.c_uint,                 # ny
                         ctypes.c_double,               # dx
                         ctypes.c_double,               # dy
                         ctypes.c_double,               # xmin
                         ctypes.c_double,               # ymin
                         ctypes.c_double,               # dt
                         ctypes.c_double,               # charge
                         ctypes.c_double,               # mass
                         ctypes.c_uint,                 # nt
                         ctypes.c_uint,                 # it0
                         ctypes.c_uint]                 # npart


        #c# xc = self._run.GetCoord(axis=0)
        #c# yc = self._run.GetCoord(axis=1)

        self._dx = self._CR['xx'][1] - self._CR['xx'][0]
        self._dy = self._CR['yy'][1] - self._CR['yy'][0]
        xc = self._CR['xx'][0] - .5*self._dx # Giving us guard cells
        yc = self._CR['yy'][0] - .5*self._dy

        func(self.r,
             self.v,
             self._E,
             self._B,
             self._B[0,:,:].shape[0],
             self._B[0,:,:].shape[1],
             self._dx,
             self._dy,
             xc, yc,
             self.dt,
             self._charge,
             self._mass,
             self._nt,
             self._it0,
             self._npart)
    #==========================================================

    #==========================================================
    def _moveall3D(self):
        """move all particles backward and forward as necessary

        Exemple  :

        Creation : 2013-05-04 15:40:32.257893

        """
        func          = self._lib.moveall3D
        func.restype  = None
        func.argtypes = [ndpointer(ctypes.c_double),    # pos
                         ndpointer(ctypes.c_double),    # vel
                         ndpointer(ctypes.c_float),     # E
                         ndpointer(ctypes.c_float),     # B
                         ctypes.c_uint,                 # nx
                         ctypes.c_uint,                 # ny
                         ctypes.c_uint,                 # nz
                         ctypes.c_double,               # dx
                         ctypes.c_double,               # dy
                         ctypes.c_double,               # dz
                         ctypes.c_double,               # xmin
                         ctypes.c_double,               # ymin
                         ctypes.c_double,               # zmin
                         ctypes.c_double,               # dt
                         ctypes.c_double,               # charge
                         ctypes.c_double,               # mass
                         ctypes.c_uint,                 # nt
                         ctypes.c_uint,                 # it0
                         ctypes.c_uint]                 # npart


        #c# xc = self._run.GetCoord(axis=0)
        #c# yc = self._run.GetCoord(axis=1)

        self._dx = self._CR['xx'][1] - self._CR['xx'][0]
        self._dy = self._CR['yy'][1] - self._CR['yy'][0]
        self._dz = self._CR['zz'][1] - self._CR['zz'][0]
        xc = self._CR['xx'][0] - .5*self._dx # Giving us guard cells
        yc = self._CR['yy'][0] - .5*self._dy
        zc = self._CR['zz'][0] - .5*self._dz
        
        print 'Going in!'
        func(self.r,
             self.v,
             self._E,
             self._B,
             self._B[0,:,:,:].shape[0],
             self._B[0,:,:,:].shape[1],
             self._B[0,:,:,:].shape[2],
             self._dx,
             self._dy,
             self._dz,
             xc, yc, zc,
             self.dt,
             self._charge,
             self._mass,
             self._nt,
             self._it0,
             self._npart)
    #==========================================================



    #==========================================================
    #==========================================================
    def _pfields(self):
        """interpolates the fields at the particle positions

        Exemple  :

        Creation : 2013-05-04 15:40:32.257893

        """
        func          = self._lib.pfields
        func.restype  = None
        func.argtypes = [ndpointer(ctypes.c_double),    # pos
                         ndpointer(ctypes.c_float),     # E
                         ndpointer(ctypes.c_float),     # B
                         ctypes.c_uint,                 # nx
                         ctypes.c_uint,                 # ny
                         ctypes.c_uint,                 # npart
                         ctypes.c_uint,                 # nt
                         ctypes.c_double,               # xmin
                         ctypes.c_double,               # ymin
                         ctypes.c_double,               # dx
                         ctypes.c_double,               # dy
                         ndpointer(ctypes.c_double),    # Ep
                         ndpointer(ctypes.c_double)]    # Bp


        self._dx = self._CR['xx'][1] - self._CR['xx'][0]
        self._dy = self._CR['yy'][1] - self._CR['yy'][0]
        xc = self._CR['xx'][0]# + self._dx# Giving us guard cells
        yc = self._CR['yy'][0]# + self._dy

        print

        func(self.r,
             self._E,
             self._B,
             self._E[0,:,:].shape[0],
             self._E[0,:,:].shape[1],
             self._npart,
             self._nt,
             xc,yc,
             self._dx,
             self._dy,
             self.Ep,
             self.Bp)

    #==========================================================

    #==========================================================
    def _pfields3D(self):
        """interpolates the fields at the particle positions

        Exemple  :

        Creation : 2013-05-04 15:40:32.257893

        """
        func          = self._lib.pfields3D
        func.restype  = None
        func.argtypes = [ndpointer(ctypes.c_double),    # pos
                         ndpointer(ctypes.c_float),     # E
                         ndpointer(ctypes.c_float),     # B
                         ctypes.c_uint,                 # nx
                         ctypes.c_uint,                 # ny
                         ctypes.c_uint,                 # nz
                         ctypes.c_uint,                 # npart
                         ctypes.c_uint,                 # nt
                         ctypes.c_double,               # xmin
                         ctypes.c_double,               # ymin
                         ctypes.c_double,               # zmin
                         ctypes.c_double,               # dx
                         ctypes.c_double,               # dy
                         ctypes.c_double,               # dz
                         ndpointer(ctypes.c_double),    # Ep
                         ndpointer(ctypes.c_double)]    # Bp


        self._dx = self._CR['xx'][1] - self._CR['xx'][0]
        self._dy = self._CR['yy'][1] - self._CR['yy'][0]
        self._dz = self._CR['zz'][1] - self._CR['zz'][0]
        xc = self._CR['xx'][0]# + self._dx# Giving us guard cells
        yc = self._CR['yy'][0]# + self._dy
        zc = self._CR['zz'][0]# + self._dy

        print

        func(self.r,
             self._E,
             self._B,
             self._E[0,:,:,:].shape[0],
             self._E[0,:,:,:].shape[1],
             self._E[0,:,:,:].shape[2],
             self._npart,
             self._nt,
             xc,yc,zc,
             self._dx,
             self._dy,
             self._dz,
             self.Ep,
             self.Bp)

    #==========================================================

    #==========================================================
    #==========================================================
    def _velinit(self, v, dv):
        """uniform random velocities

        Creation : 2013-04-28 17:12:42.773824

        """

        # create the array
        self._v0 = np.zeros((3, self._npart))

        # x, y, z components of the velocity
        c = 0

        # loop over vx0, vy0, vz0
        for v_i,dv_i in zip(v,dv):

            v1 = v_i - dv_i/2.
            v2 = v_i + dv_i/2.

            self._v0[c,:] = np.random.random(self._npart)
            self._v0[c,:] = self._v0[c,:] * (v2 - v1) + v1
            c += 1

    #==========================================================

    #==========================================================
    #==========================================================
    def _velinitn(self, v, dv):
        """uniform random velocities

        Creation : 2015-10-09 17:12:42.773824

        """

        # create the array
        self._v0 = np.zeros((3, self._npart))

        # x, y, z components of the velocity
        c = 0

        # loop over vx0, vy0, vz0
        for v_i,dv_i in zip(v,dv):

            self._v0[c,:] = np.random.randn(self._npart)*dv_i
            self._v0[c,:] = self._v0[c,:] + v_i
            c += 1

    #==========================================================
    #==========================================================
    def _velinit_mag(self, v, dv):
        """uniform random velocity magnetudies but constant unit vector

        Creation : 2015-04-28 17:12:42.773824

        """

        # create the array
        self._v0 = np.zeros((3, self._npart))

        # x, y, z components of the velocity
        c = 0

        # loop over vx0, vy0, vz0
        vmag = 0.0
        for v_i in v:
            vmag = vmag+v_i**2
        vmag = sqrt(vmag)
        rmag = np.random.random(self._npart)
        for v_i in v:
            self._v0[c,:] = dv*rmag*v_i/vmag
            c += 1

    #==========================================================
    #==========================================================
    def load_randu(self,
                   r,
                   dr,
                   v,
                   dv):

        """ loads particles randomly in a rectangle

        The method loads the particles in a rectangle
        of side 'dr' entered around 'r'.
        Particles are loaded uniformly in that rectangle.
        Velocities are uniformly distributed within 'dv' from 'v'

        Carefull : this method will erase any initial position
        and velocity that may have already been initialized.

        @param r    : (x,y,z) center of the rectangle
        @param dr   : (dx,dy,dz) size of the rectangle
        @param v    : (Vx,Vy,Vz) mean initial velocity
        @param dv   : (dVx,dVy,dVz) initial velocity interval

        @return: @todo

        Exemple  :

        Creation : 2013-05-01 11:21:19.671182

        """
        print 'Loading particles in a spatially uniform random distribution...'

        self._r0  = np.zeros((3, self._npart))

        c = 0. # x, y and z components

        for r_i, dr_i in zip(r, dr):

            r1            = r_i - dr_i/2.
            r2            = r_i + dr_i/2.
            self._r0[c,:] = np.random.random(self._npart)

            self._r0[c,:]  = self._r0[c,:]* (r2 - r1) + r1

            c += 1

        self._velinit(v,dv)
        #self._velinitn(v,dv)

    #==========================================================

    #==========================================================
    #==========================================================
    def load_randp(self,
                   r,
                   dr,
                   v,
                   dv):

        """ loads particles randomly in a rectangle

        The method loads the particles in a rectangle
        of side 'dr' entered around 'r'.
        Particles are loaded uniformly in that rectangle.
        Velocities are uniformly distributed within 'dv' from 'v'

        Carefull : this method will erase any initial position
        and velocity that may have already been initialized.

        @param r    : (x,y,z) center of the rectangle
        @param dr   : (dx,dy,dz) size of the rectangle
        @param v    : (Vx,Vy,Vz) mean initial velocity
        @param dv   : (dVx,dVy,dVz) initial velocity interval

        @return: @todo

        Exemple  :

        Creation : 2013-05-01 11:21:19.671182

        """
        print 'Loading particles in a spatially uniform random distribution...'

        self._r0  = np.zeros((3, self._npart))

        c = 0. # x, y and z components

        for r_i, dr_i in zip(r, dr):

            r1            = r_i - dr_i/2.
            r2            = r_i + dr_i/2.
            self._r0[c,:] = np.random.random(self._npart)

            self._r0[c,:] = self._r0[c,:]* (r2 - r1) + r1

            c += 1

        #self._velinit(v,dv)
        self._velinitn(v,dv)

    #==========================================================

    #==========================================================
    #==========================================================
    def load_randu_mag(self,
                      r,
                      dr,
                      v,
                      dv):

        """ loads particles randomly in a rectangle

        The method loads the particles in a rectangle
        of side 'dr' entered around 'r'.
        Particles are loaded uniformly in that rectangle.
        Velocities are uniformly distributed within 'dv' from 'v'
        However dv is now broken into two parts, dv|| and dv perp

        Carefull : this method will erase any initial position
        and velocity that may have already been initialized.

        @param r    : (x,y,z) center of the rectangle
        @param dr   : (dx,dy,dz) size of the rectangle
        @param v    : (Vx,Vy,Vz) mean initial velocity
        @param dv   : (dV||,dVperp) initial velocity interval

        @return: @todo

        Exemple  :

        Creation : 2013-05-01 11:21:19.671182

        """
        print 'Loading particles in a spatially uniform random distribution...'

        self._r0  = np.zeros((3, self._npart))

        c = 0. # x, y and z components

        for r_i, dr_i in zip(r, dr):

            r1            = r_i - dr_i/2.
            r2            = r_i + dr_i/2.
            self._r0[c,:] = np.random.random(self._npart)

            self._r0[c,:]  = self._r0[c,:]* (r2 - r1) + r1

            c += 1

        self._velinit_mag(v,dv)

    #==========================================================





    #==========================================================
    #==========================================================
    def load_randn(self,
                   r,
                   dr,
                   v,
                   dv):

        """ loads particles randomly in a spatial gaussian

        The method loads the particles in a spatial gaussian
        of std 'dr' entered around 'r'.
        Particles are loaded in a normal distribution.
        Velocities are uniformly distributed within 'dv' from 'v'

        Carefull : this method will erase any initial position
        and velocity that may have already been initialized.

        @param r  : (x,y,z) center of the gaussian
        @param dr : (dx,dy,dz) standard deviation from r
        @param v  : (Vx,Vy,Vz) mean initial velocity
        @param dv : (dVx,dVy,dVz) initial velocity interval

        @return: @todo

        Exemple  :

        Creation : 2013-05-01 11:21:19.671182

        """
        print 'Loading particles in a spatially gaussian random distribution...'

        self._r0  = np.zeros((3, self._npart))

        c = 0. # x, y and z components

        for r_i, dr_i in zip(r, dr):

            r1            = r_i - dr_i/2.
            r2            = r_i + dr_i/2.
            self._r0[c,:] = np.random.random(self._npart)

            self._r0[c,:]  = self._r0[c,:]* (r2 - r1) + r1

            c += 1

        self._velinit(v,dv)
