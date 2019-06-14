#######################################################################
#                                                                     #
#                  Python Progs :  vdist_plotter.py                   #
#                  Aruthor      :  Colby Haggerty                     #
#                  Date         :  2016.11.13                         #
#                                                                     #
#                                                                     #
#######################################################################
import datetime
import numpy as np
import matplotlib.pylab as plt
import pdb
from scipy.ndimage import gaussian_filter
from .dumpID import DumpID
from .vdist import VDist

class VDistPlotter(object):
    """ velocity distribution fucntion plotter
    """
#===========================================================#

    def __init__(self,
                 r0=None,
                 dx=None,
                 **dumpargs): 
        """ Initilazition Routine for the p3d_run object

            **dumpargs are the kwargs for the DumpID object
                 They include and are limited to the args:
                 num, param_file and path
        """

        self.update_box(r0, dx)
        self._dumpid = DumpID(**dumpargs)

#===========================================================#

    def update_box(self, r0, dx):
        self.r0    = self._get_r0(r0)
        self.dx    = self._get_dx(dx)
        self.parts = None
        self._par  = None
        return self

#===========================================================#

    def plot2d(self,
               d1,
               d2, 
               sp,
               ax=None,
               dz=None,
               par=False,
               v0_frame=False,
               smooth=0.,
               ctargs=None,
               pcmargs={},
               **kwargs):

        """ Plots a 2D distro
        Args:
            d1 (int 0-2): The velocity space direction coresponding to the x
                axis

            d2 (int 0-2): The velocity space direction coresponding to the y
                axis
                Note that if par == False:
                    0 -> x, 1 -> y, 2 -> z
                And if par == True:
                    0 -> b, 1 -> exb, 2 -> bx(exb)

            ax (Matplotlib axes obj, optional): the subplot you want this to 
                be plotted to. If left as None, the code grabs the current axes

            sp (str 'i' or 'e'): the species that you want plotted

            dz (float, optional): The width of integrated "pizza" in the 3rd 
                velocity space direction. Default is None (All space)

            v0_frame (boo, optional): If true, this shifts the center of the
                pizza to the mean velocity in the 3rd direction. Default is 
                set to False

            smooth (float, optional) the std for the guassian filter to smooth
                the histogram

            ctargs (dict): Keyword arguments for the contour function

            pcmargs (dict): Keyword arguments for the pcolormesh function

            **kwargs: Keyword arguments that will ultimtily be passed to the
                numpy histogram2d function, so look there for more details

        Returns:
            ax: Retuns matplotlib AxesSubplot of coresponding to plotted fig
        """
        from matplotlib.colors import LogNorm
        
        ax = self._get_ax(ax)
        self._set_parts(par)
        # Debug stuff
        k1,k2,k3 = self._get_coord(d1,d2)
        v1 = self.parts[sp][k1]
        v2 = self.parts[sp][k2]
        if dz is None:
            v3 = None
        else:
            v3 = self.parts[sp][k3]

        H, xx, yy = VDist().vdist2d(v1, v2, v3, dz, v0_frame, **kwargs)
        H = self._smooth(xx, yy, H, smooth)

        pcmargs = dict(pcmargs)
        
        if pcmargs.pop('norm',None):
            pcmargs['norm'] = LogNorm(vmax=pcmargs.pop('vmax',None),
                                      vmin=pcmargs.pop('vmin',None))
        #pdb.set_trace() 
        pcm = ax.pcolormesh(xx, yy, H, **pcmargs)
        ax.set_aspect('equal')

        if ctargs is not None:
            self._add_contours(ax, xx, yy, H, ctargs)

        self._set_labels(ax, k1, k2, k3, sp, dz)
        self._set_lims(ax, xx, yy, pcm)

# For some reason this is just not working right now!
        self._set_minorticks_on(ax)
        #pdb.set_trace()

        return ax, pcm

#===========================================================#

    def spec1d(self,
               dir,
               pitch_angle,
               delta_pitch,
               species,
               mass=None,
               ax=None,
               v0_frame=False,
               v_light=None,
               smooth=0.,
               ctargs=None,
               pcmargs={},
               **kwargs):

        """ Plots a 2D distro
        Args:
            pitch angle (int 0-2): The velocity space direction coresponding to the x
                axis
            dpitch_angle (Matplotlib axes obj, optional): the subplot you want this to 
                be plotted to. If left as None, the code grabs the current axes
            species (str 'i' or 'e'): the species that you want plotted
            dz (float, optional): The width of integrated "pizza" in the 3rd 
                velocity space direction. Default is None (All space)
            v0_frame (bool, optional): If true, this shifts the center of the
                pizza to the mean velocity in the 3rd direction. Default is 
                set to False
            ctargs (dict): Keyword arguments for the contour function
            pcmargs (dict): Keyword arguments for the pcolormesh function
            **kwargs: Keyword arguments that will ultimtily be passed to the
                numpy histogram2d function, so look there for more details

        Returns:
            ax: Retuns matplotlib AxesSubplot of coresponding to plotted fig
        """
        
        ax = self._get_ax(ax)
        self._set_parts(par=True)

        if dir in [0,1,2]:
            dir = 'xyz'[dir]

        if species == 'i':
            mass = 1.
        elif species == 'e':
            if mass is None:
                mass = float(raw_input('Enter electron mass: '))


        H, xx, yy = VDist().spec1d(self.parts[species],
                                   dir=dir,
                                   pa=pitch_angle, 
                                   dpa=delta_pitch,
                                   mass=mass,
                                   v0_frame=v0_frame,
                                   v_light=v_light,
                                   **kwargs)

        #H = self._smooth(xx, yy, H, smooth)

        #pdb.set_trace() 
        pcm = ax.pcolormesh(xx, yy, H, **pcmargs)
        #ax.set_aspect('equal')

        if ctargs is not None:
            raise NotImplementedError()
            #self._add_contours(ax, xx, yy, H, ctargs)

        #self._set_labels(ax, k1, k2, k3, sp, dz)
        self._set_lims(ax, xx, yy, pcm)
        ax.set_xlabel(dir+' $(d_i)$')
        ax.set_ylabel('KE $(m_iv_0^2)$')
        ax.set_title('$\\theta = [{}, {}]$'.format(pitch_angle-delta_pitch/2., pitch_angle+delta_pitch/2.))

# For some reason this is just not working right now!
        self._set_minorticks_on(ax)
        #pdb.set_trace()

        return ax, pcm
        
#===========================================================#

    def _add_contours(self, ax, xx, yy, H, ctargs):
        if type(ctargs) is not dict:
            msg = 'ctargs is type {}, but must be type dict.'\
                  'No conrours plotted!'
            print msg.format(type(ctargs))
            pass 
        
        defaults = dict(linestyles='solid', colors='black')
        for k in defaults:
            ctargs[k] = ctargs.get(k,defaults[k])
    
        redge = lambda x : (x[1:] + x[:-1])/2.
        
        ctr = ax.contour(redge(xx) , redge(yy), H, **ctargs)

        return ctr


#===========================================================#

    def _smooth(self, xx, yy, H, smooth):
        if smooth:
            from scipy.ndimage import gaussian_filter as gf

            # smooth is assumed to be in velocity units
            # So we need to convert it to grid units
            # Also note that we are assuming dvx = dvy which does 
            # not have to be true!
            sig = smooth/(xx[1] - xx[0] + yy[1] - yy[0])*2.

            #return gf(H, sigma=sig, mode='constant', cval=0.)
            return gf(H, sigma=sig)

        else:
            return H

#===========================================================#

    def _set_minorticks_on(self, ax):
        temp_ax = plt.gca()
        plt.sca(ax)
        plt.minorticks_on()
        plt.sca(temp_ax)

        return self

#===========================================================#

    def _set_lims(self, ax, xx, yy, pcm):
        from matplotlib.colors import LogNorm

        ax.set_xlim(xx[0], xx[-1])
        ax.set_ylim(yy[0], yy[-1])
        #pcmax = pcm.get_clim()[1]
        #pcm.set_clim(pcmax/10., pcmax)
        #pcm.set_norm(LogNorm())

#===========================================================#

    def _set_labels(self, ax, k1, k2, k3, sp, dz):
        label_dict = {'vx': '$V_x$',
                      'vy': '$V_y$',
                      'vz': '$V_z$',
                      'v0': '$V_\parallel$',
                      'v1': '$V_{E\\times B}$',
                      'v2': '$V_{B\\times(E\\times B)}$'}
        
        ax.set_xlabel(label_dict[k1])
        ax.set_ylabel(label_dict[k2])

        title = '"{}", '.format(sp) 
        for v,r0,dx in zip('XYZ', self.r0, self.dx):
            title += '{}[{:.2f}, {:.2f}], '.format(v, r0-dx/2., r0+dx/2.)
            if v == 'Y':
                title += '\n'

        if dz is not None:
            title += ', d{} [{:.2f}]'.format(label_dict[k3], dz)

        ax.set_title(title,size=8)
        return self

#===========================================================#

    def _get_coord(self, d1, d2):
        if d1 not in range(3) and d2 not in range(3):
            msg = '{0} or {1} not possible option (0, 1, 2).\n'\
                  'Remember (0, 1, 2) -> (x, y, z) or \n'\
                  'if par == True (0, 1, 2) -> (b, exb, bx(exb))'
            print msg.format(d1,d2)
            raise KeyError('Bad dimension numbers {} or {}'.format(d1,d2))
        elif d1 == d2:
            msg = '{0} = {1} which seems like a mistake.\n'\
                  'Also the 3rd direction is ambiguous\n'\
                  'Remember (0, 1, 2) -> (x, y, z) or \n'\
                  'if par == True (0, 1, 2) -> (b, exb, bx(exb))'
            print msg.format(d1,d2)
            raise KeyError('Bad dimension numbers: {} = {}'.format(d1,d2))

        if self._par:
            coord_map = {0: 'v0',
                         1: 'v1',
                         2: 'v2'}
        else:
            coord_map = {0: 'vx',
                         1: 'vy',
                         2: 'vz'}
       
        d3 = range(3)
        d3.remove(d1)
        d3.remove(d2)
        d3 = d3[0]
        
        return coord_map[d1], coord_map[d2], coord_map[d3]
        
    
#===========================================================#

    def _set_parts(self, par):
        if self.parts is None or self._par != par:
            self._par = par
            self.parts = self._dumpid.get_part_in_box(
                self.r0, self.dx, self._par)

        return self

#===========================================================#

    def _get_r0(self, r0):
        atmp_tol = 5
        ctr  = 0
        msg = 'Please enter the center location to \n'\
              'calculate the distribution function (x, y, (z)): '
        while not isinstance(r0, (list, np.ndarray, np.generic)):
            r0 = raw_input(msg).split(',')
            r0 = self._convert_to_float(r0)

            if ctr > atmp_tol:
                raise IOerror('Exceed allowed attemps')
            else:
                ctr += 1

        return r0

#===========================================================#

    def _get_dx(self, dx):
        atmp_tol = 5
        ctr  = 0
        msg = 'Please enter the width of the box to collect \n'\
               'particles in  the distribution function (dx, dy, (dz)): '
        while not isinstance(dx, (list, np.ndarray, np.generic)):
            dx = raw_input(msg).split(',')
            dx = self._convert_to_float(dx)

            if ctr > atmp_tol:
                raise IOerror('Exceed allowed attemps')
            else:
                ctr += 1

        return dx

#===========================================================#

    def _get_ax(self, ax):
        if ax is None:
            #print "Grabbing matplotlib's current axes"
            ax = plt.gca()

        return ax

#===========================================================#

    def _convert_to_float(self, ar):
        br = np.zeros(len(ar))
        for c,r in enumerate(ar):
            try:
                br[c] = float(r)
            except ValueError:
                return None

        return br
            
                    
