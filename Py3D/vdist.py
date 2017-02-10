#######################################################################
#                                                                     #
#                  Python Progs :  vdist.py                           #
#                  Aruthor      :  Colby Haggerty                     #
#                  Date         :  2016.02.08                         #
#                                                                     #
#                                                                     #
#######################################################################
import os
import sys 
import datetime
import numpy as np
import struct
import glob
import pdb
import warnings
from scipy.io.idl import readsav
from scipy.ndimage import gaussian_filter

class VDist(object):
    """ velocity distribution fucntion calculator
    """

    def __init__(self): 
        """ Initilazition Routine for the p3d_run object
        """

        pass

    def vdist2d(self,
                v1,
                v2,
                v3=None,
                dz=None,
                v0_frame=False,
                v_light=None,
                **kwargs):

        """ Simple 2D velocity histogram
        """
        
        if v3 is not None:
            v1, v2 = self._trim_v3(v1, v2, v3, dz, v0_frame)
        
        if kwargs.get('bins') is None and kwargs.get('range') is None:
           self._set_bins(v1,v2, **kwargs)

        H, x_edge, y_edge = np.histogram2d(v1, v2, **kwargs)

        return H.T, x_edge, y_edge
        

    def _trim_v3(self, v1, v2, v3, dz, v0_frame):

        if dz is None:
            dz = np.std(v3)
        if v0_frame:
            v3 = v3 - np.mean(v3)
        
        ind = np.where(abs(v3) < dz/2.)
        v1 = v1[ind]
        v2 = v2[ind]

        return v1, v2


    def _set_bins(self, v1, v2, **kwargs):

        kwargs.pop('range', None)
        kwargs.pop('bins', None)

        v1_mean = np.mean(v1)
        v2_mean = np.mean(v2)

        v_std = np.std(np.sqrt(v1**2 + v2**2))
        v_nsteps = 50
        # dv = v_std/10.

        kwargs['bins'] = [np.linspace(v1_mean - 2.5*v_std,
                                      v1_mean + 2.5*v_std,
                                      v_nsteps),
                          np.linspace(v2_mean - 2.5*v_std,
                                      v2_mean + 2.5*v_std,
                                      v_nsteps)]



    def vdist2d_pitch(self,
                      v1,
                      v2,
                      v3,
                      pa=90.,
                      dpa=20.,
                      v0_frame=False,
                      **kwargs): 
        """ Pitch party!!!!
            
            pa = pitch angle
            dpa = angular width around pa (Delta Pitch Angle)
        """

        if v0_frame:
            v1 = v1 - np.mean(v1)
            v2 = v2 - np.mean(v2)
            v3 = v3 - np.mean(v3)

        pitch_angle = -1.0*(np.arctan(v3/np.sqrt(v1**2+v2**2))/np.pi*180. - 90.)
        ind = np.where(abs(pitch_angle - pa) < dpa/2.)
        
        v1 = v1[ind]
        v2 = v2[ind]
        
        if kwargs.get('bins') is None and kwargs.get('range') is None:
           self._set_bins(v1, v2, **kwargs)

        H, x_edge, y_edge = np.histogram2d(v1, v2, **kwargs)

        H = H.T

        xx,yy = np.meshgrid(x_edge, y_edge)

        norm_cone = self._int_cone(pa, dpa, x_edge, y_edge)

        H = H/norm_cone

        return H, x_edge, y_edge


    def _get_gamma(self, v1, v2, v3, v_light):
        v2 = v1**2 + v2**2 + v3**2
        c2 = v_light**2


 def spec1d(pts,
            dir,
            pa,
            dpa,
            mass,
            v_light=None,
            v0_frame=False,
            **kwargs):
    """ core program for calculating the energy spectrum along a direction"""

    if v_light is not None:
        gamma = self._get_gamma(pts['vx'], pts['vy'], pts['vz'], v_light)
        c2 = v_light**2

        KE = mass*c2*(gamma - 1.)
    else:
        KE = .5*mass*(pts['vx']**2 + pts['vy']**2 + pts['vz']**2)
    
    v0 = pts['v0']
    v1 = pts['v1']
    v2 = pts['v2']
    if v0_frame:
        v0 = v0 - np.mean(v0)
        v1 = v1 - np.mean(v1)
        v2 = v2 - np.mean(v2)

    # Note: Does pitch angle change in special relativity? 
    #       This might be pretty wrong. Maybe ask someone?
    pitch_angle = -1.0*(np.arctan(v0/np.sqrt(v1**2 + \
        v2**2))/np.pi*180. - 90.)

    ind = np.where(abs(pitch_angle - pa) < dpa/2.)

    H, x_edge, y_edge = np.histogram2d(pts[dir][ind],
                                       KE[ind],
                                       normed=True,
                                       **kwargs)

    eng = (y_edge[:-1] + y_edge[1:])/2.

    if v_light is not None:
# I stole this from EFlux code, I dont remeber why we do this
# so you know, use at your own risk
        warn_msg = 'Relativistic energy code implemented without though:' +
                   'Use with caution!!!'
        warnings.warn(msg)
        Ebar = eng/mass/v_light**2
        rel_vel = np.sqrt((Ebar**2 + 2.*Ebar)/(Ebar**2 + 2.*Ebar + 1.))
        rel_vel = rel_vel*v_light
        H = (H*eng*rel_vel).T

    else:
        H = (H*eng**2/sqrt(eng)).T


    return H, x_edge, y_edge

        

    def eflux(self, v1, v2, v3, mass, v_light=None, **kwargs):
        """ Function that makes an energy flux distrobution
            ToDo: add an eflux_autobin funciton
        """

        pitch_angle = -1.0*(np.arctan(v3/np.sqrt(v1**2+v2**2))/np.pi*180. - 90.)

        v2 = v1**2 + v2**2 + v3**2

        if v_light is not None:
            gamma = self._get_gamma(v1, v2, v3, v_light)
            c2 = v_light**2

            KE = mass*c2*(gamma - 1.)
        else:
            KE = .5*mass*v2

        H, x_edge, y_edge = np.histogram2d(KE, pitch_angle, normed=True, **kwargs)
        H = H.T

        xx,yy = np.meshgrid(x_edge,y_edge)
        eng = (xx[1:,1:] + xx[:-1,:-1])/2.

        cos_norm = (np.cos(yy[:-1,:-1]/180.*np.pi) - np.cos(yy[1:,1:]/180.*np.pi))
        cos_norm = cos_norm/sum(cos_norm)

# We need to account for a geometric factor because
# E is basicly v^2, so we are binning circles?

        
        if v_light is None:
            H = H/np.sqrt(eng)
            H = H/cos_norm*eng**2.*np.size(pitch_angle)

        else:
            Ebar = eng/mass/v_light**2
            rel_vel = np.sqrt((Ebar**2 + 2.*Ebar)/(Ebar**2 + 2.*Ebar + 1.))
            rel_vel = rel_vel*v_light

            H = H/cos_norm*eng*rel_vel*np.size(pitch_angle)

        return H, x_edge, y_edge

    def _int_cone(self,pa,dpa,xedges,yedges):
        """
#-----------------------------------------------------------------------------
#   Method      : _int_cone
#
#   Discription : This integrates a cone over a differental grid
#
#   Args        : xedges the x edges of the histogram
#               : yedges the y edges of the histogram
#
#   Comments    : I think this is working ok? 
#-----------------------------------------------------------------------------
        """
        

        xx,yy = np.meshgrid(yedges,xedges)
        intcone = 1./18.*(
          (-xx[:-1,:-1]**3+6.*xx[:-1,:-1]*yy[:-1,:-1]*np.sqrt(xx[:-1,:-1]**2+
           yy[:-1,:-1]**2) + 3.*yy[:-1,:-1]**3*np.log(np.sqrt(xx[:-1,:-1]**2+
           yy[:-1,:-1]**2) + xx[:-1,:-1] + np.spacing(4)) +
           3.*xx[:-1,:-1]**3*np.log(np.sqrt(xx[:-1,:-1]**2+
           yy[:-1,:-1]**2)+yy[:-1,:-1] + np.spacing(4))) -
          (-xx[1:,1:]**3+6.*xx[1:,1:]*yy[:-1,:-1]*np.sqrt(xx[1:,1:]**2+
           yy[:-1,:-1]**2) + 3.*yy[:-1,:-1]**3*np.log(np.sqrt(xx[1:,1:]**2+
           yy[:-1,:-1]**2) + xx[1:,1:] + np.spacing(4)) + 
           3.*xx[1:,1:]**3*np.log(np.sqrt(xx[1:,1:]**2+
           yy[:-1,:-1]**2)+yy[:-1,:-1] + np.spacing(4))) -
          (-xx[:-1,:-1]**3+6.*xx[:-1,:-1]*yy[1:,1:]*np.sqrt(xx[:-1,:-1]**2+
           yy[1:,1:]**2) + 3.*yy[1:,1:]**3*np.log(np.sqrt(xx[:-1,:-1]**2+
           yy[1:,1:]**2) + xx[:-1,:-1] + np.spacing(4)) + 
           3.*xx[:-1,:-1]**3*np.log(np.sqrt(xx[:-1,:-1]**2+
           yy[1:,1:]**2)+yy[1:,1:] + np.spacing(4))) +
          (-xx[1:,1:]**3+6.*xx[1:,1:]*yy[1:,1:]*np.sqrt(xx[1:,1:]**2+
           yy[1:,1:]**2) + 3.*yy[1:,1:]**3*np.log(np.sqrt(xx[1:,1:]**2+
           yy[1:,1:]**2) + xx[1:,1:] + np.spacing(4)) + 
           3.*xx[1:,1:]**3*np.log(np.sqrt(xx[1:,1:]**2+
           yy[1:,1:]**2)+yy[1:,1:] + np.spacing(4))))

        norm = 1./18.*(
               (-xx[0,0]**3+6.*xx[0,0]*yy[0,0]*np.sqrt(xx[0,0]**2+yy[0,0]**2) + 
                3.*yy[0,0]**3*np.log(np.sqrt(xx[0,0]**2+yy[0,0]**2) + xx[0,0] + 
                np.spacing(4)) + 
                3.*xx[0,0]**3*np.log(np.sqrt(xx[0,0]**2+yy[0,0]**2)+yy[0,0] + 
                np.spacing(4))) -
               (-xx[-1,-1]**3+6.*xx[-1,-1]*yy[0,0]*np.sqrt(xx[-1,-1]**2+yy[0,0]**2) + 
                3.*yy[0,0]**3*np.log(np.sqrt(xx[-1,-1]**2+yy[0,0]**2) + xx[-1,-1] + 
                np.spacing(4)) + 
                3.*xx[-1,-1]**3*np.log(np.sqrt(xx[-1,-1]**2+yy[0,0]**2)+yy[0,0] + 
                np.spacing(4))) -
               (-xx[0,0]**3+6.*xx[0,0]*yy[-1,-1]*np.sqrt(xx[0,0]**2+yy[-1,-1]**2) + 
                3.*yy[-1,-1]**3*np.log(np.sqrt(xx[0,0]**2+yy[-1,-1]**2) + xx[0,0] + 
                np.spacing(4)) + 
                3.*xx[0,0]**3*np.log(np.sqrt(xx[0,0]**2+yy[-1,-1]**2)+yy[-1,-1] + 
                np.spacing(4))) +
               (-xx[-1,-1]**3+6.*xx[-1,-1]*yy[-1,-1]*np.sqrt(xx[-1,-1]**2+
                yy[-1,-1]**2) + 
                3.*yy[-1,-1]**3*np.log(np.sqrt(xx[-1,-1]**2+yy[-1,-1]**2) +
                xx[-1,-1] + np.spacing(4)) + 
                3.*xx[-1,-1]**3*np.log(np.sqrt(xx[-1,-1]**2+yy[-1,-1]**2)+
                yy[-1,-1] + np.spacing(4))))
        
        self.normie=norm

        print 'norm = ',norm
        print 'otha = ',abs(1./np.tan((pa-dpa/2.)/180.*np.pi) - 1./np.tan((pa+dpa/2.)/180.*np.pi))

        #return intcone/norm#*abs(1./np.tan((pa-dpa/2.)/180.*np.pi) - 1./np.tan((pa+dpa/2.)/180.*np.pi))
        return intcone/intcone.min()#*abs(1./np.tan((pa-dpa/2.)/180.*np.pi) - 1./np.tan((pa+dpa/2.)/180.*np.pi))


