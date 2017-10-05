from scipy.ndimage import gaussian_filter as gf
import numpy as np
from .movie import Movie
from .sub import ims
from .sub import rotate_ten
from .sub import calc_psi
from .sub import date_file_prefix
import pdb
#import _methods

class PatPlotter(object):
    """ A class to load and plot a simulation in the style of kittypat
    """

    def __init__(self, name_style='p3d', **mvargs):
        """ Make a set of plots like kittypat used to

        :param mvargs: Keyword arguments inteneded to be passed to py3d.Movie
        :type mvargs: kwargs or dict
        """

        self._M = Movie(**mvargs)

        self.ctrs = []

        pgv2 = ['ni rho bx ex by ey bz ez |b| jx jy jz',
               'jix jex jiy jey jiz jez vix vex viy vey viz vez']
        pgv1 = ['bx by bz |b|', 'ni ne rho',
                'ni ne', 'piyy peyy',
                'ex ey ez', 'jx jy jz',
                'vix viy viz', 'vex vey vez',
                'tixx tiyy tizz', 'texx teyy tezz',
                'tipar tiperp1', 'tepar teperp1']

        self.page_vars_2D = [k.split() for k in pgv2]
        self.page_vars_1D = [k.split() for k in pgv1]

#========================================

    def make_plots(self, 
                   time=None,
                   slc=None,
                   xy_lims=None,
                   cut_dir='y',
                   cut_locs=None,
                   cont_override=None,
                   **kwargs):

        if slc is not None:
            slc = slc[::-1] 

        # These will need to be arguments
        self.sig = 0
        #self.xy_lims = [[4.5, 6.], [7.3, 8.]]

        # Commented out
        self.d  = self._load_movie(time, slc)
        self._cut_dir =  cut_dir
        self._not_cut_dir = 'xy'[abs(1 - 'xy'.find(self._cut_dir))]
        self._created_files = []

        # This is gonna need a whole thing
        psi0, levels = self._calc_psi0(cont_override)
        self.psi0 = psi0
        
        ctargs = dict(linewidths=.5, levels=levels)

        #Skipping this too for right now
        self._ctrs = self._gen_conts(**ctargs)

        ldgargs = dict(loc=3, ncol=10, borderaxespad=0., frameon=0, 
            prop={'size':6}, bbox_to_anchor=(0., 1.02, 1., .102))
        
        # This needs to be automated
        _x,_i = self._gen_xcuts(cut_locs)
        self.xpcs = _x
        self.ipcs = _i

        self.page_counter = 0
        for page_vars in self.page_vars_2D:
            self.clean_fig()
            self._plot2D(page_vars, **kwargs)
            self._plot2D_cutlines()
            self._plot2D_set_lims(xy_lims)
            self.savefig()

        for ipc in self.ipcs:
            self.clean_fig()
            self._plot1D(ipc, ldgargs=ldgargs)
            self.plot_psi_intercepts(ipc)
                
            self.savefig()

        self._make_pdf()

#========================================

    def _plot2D(self, page_vars, var_labels=None, **kwargs):

        pcm = []

        if var_labels is None: var_labels = page_vars

        for a,v,l in zip(self.ax, page_vars, var_labels):

            print 'Plotting {}...'.format(v)
            vrs = gf(self.d[v], sigma=self.sig)
            pcm += [ims(self.d, vrs, ax=a, no_draw=1, **kwargs)]

            for c in self._ctrs:
                a.plot(*c, color='k',linewidth=.5)
            
            _title = l + ':{:.3f} {:.3f}'.format(vrs.min(), vrs.max())
            a.set_title(_title, fontsize=6)

            if a not in self.fig.axes[-2:]:
                a.set_xlabel('')

        return pcm

#========================================

    def _plot2D_cutlines(self):
        lines = []
        for a in self.ax:
            for ipc in self.ipcs:
                lines += a.plot(self.d['yy']*0. + self.d['xx'][ipc], 
                                self.d['yy'],
                                'k--',
                                linewidth=.5)
        return lines

#========================================

    def _plot2D_set_lims(self, xy_lims):
        if xy_lims:
            for a in self.ax:
                a.set_xlim(xy_lim[0])
                a.set_ylim(xy_lim[1])


#========================================

    def _gen_xcuts(self, cut_locs):
        ll = (self.d[2*self._not_cut_dir][-1] +
              self.d[2*self._not_cut_dir][0])
        
        if cut_locs is None:
            ncuts = 10
            #xpcs = np.arange(4.6, 6., .1)
            xpcs = np.arange(10)/10.*ll

        else:
            xpcs = cut_locs

        ipcs = [abs(self.d['xx'] - k).argmin() for k in xpcs]

        return xpcs, ipcs

#========================================

    def _gen_conts(self, **ctargs): # d should containe psi
        d = self.d
        print 'Generating Contours...'
        import matplotlib.pyplot as plt
        plt.ioff()
        _f,a = plt.subplots(1,1,num=1)

        ctargs['colors'] = ctargs.get('colors','k')
        ctargs['linestyles'] = ctargs.get('linestyles', 'solid')

        cts = a.contour(d['xx'], d['yy'], d['psi'].T, **ctargs)

        ctrs = []
        for c in cts.collections:
            for k in c.get_paths():
                ctrs.append(k.vertices.T)
        
        _f.clf()
        del(_f)

        return ctrs

#========================================

    def _plot1D(self, ip, var_labels=None,
               ptargs=None,
               ldgargs=None,
               xlim=None,
               ylim=None):
        """ Function to make a plot of a 1D cut

            Parameters
            =========
            ax : (matplotlib.pyplot.axis)
            xy : (numpy.array)
            vrs : [numpy.array]

            Todo: Fix the stuipd kwargs
        """

        if var_labels is None: var_labels = self.page_vars_1D

        # We will come back it this if we have to
        #if not ptargs: 
        #    ptargs = 3*({},)
        #elif type(ptargs) is dict:
        #    ptargs =3*(ptargs,)

        _xy = self.d[2*self._cut_dir]
        _lim = _xy[[0,-1]]
        _cut_index = np.s_[ip,:]
        if self._cut_dir is 'x': _cut_index = _cut_index[::-1]

        lines = []
        for a,vrs,labs in zip(self.ax, self.page_vars_1D, var_labels):

            _sl = []
            #for v,l,pwargs in zip(vrs, var_labels, plot_kwargs):
            for v,l in zip(vrs, labs):
                gfv = gf(self.d[v][_cut_index], sigma=self.sig)
                _sl += a.plot(_xy, gfv, label=l)
        
            lines += [_sl]

            a.set_xlim(_lim)
            
            _tl = 'cut @ {} = {:1.2f}'.format(self._not_cut_dir, 
                self.d[2*self._not_cut_dir][ip])

            a.set_title(_tl, size=6, loc='right')

            if ldgargs:
                if type(ldgargs) is not dict: ldgargs = {}
                a.legend(**ldgargs)
            
            a.minorticks_on()

        return lines

#========================================

    def _calc_midplane(self):
        return np.abs(self.d['bx']).argmin(axis=1)

#========================================

    def _calc_psi0(self, cont_override):
        if cont_override is None:
            mp = self._calc_midplane()
            xrng = np.arange(len(mp))

            if self.d['bx'][0,-1] - self.d['bx'][0,0] > 0.:
                arg_psi0 = self.d['psi'][xrng, mp].argmax()
            else:
                arg_psi0 = self.d['psi'][xrng, mp].argmin()
            
            psi0 = self.d['psi'][xrng[arg_psi0], mp[arg_psi0]]
            npsi = 10
            dpsi = np.abs(psi0 - self.d['psi'][arg_psi0, 4])/(npsi/2.)

            levels= np.arange(psi0 - dpsi*npsi, psi0 + dpsi*npsi, dpsi)
       
        else:
            ctovrd_msg = ('If cont_override is not None then it '
                          'must be (psi0, levels)\n'
                          'i.e. tuple(float, list (or numpy.array))')
            assert type(cont_override) is tuple, ctovrd_msg
            assert type(cont_override[0]) in [float, np.float64], ctovrd_msg
            assert type(cont_override[1]) in [list,  np.ndarray], ctovrd_msg

            psi0, levels = cont_override

        return psi0, levels

#========================================
    def plot_psi_intercepts(self, ipc):
        
        for a in self.ax:
            yl = a.get_ylim()
            for kp in self._psi_intercepts(ipc):
                xl = 2*[self.d[2*self._cut_dir][kp]]
                a.plot(xl, yl, 'k--', linewidth=.5)

            a.set_ylim(yl)

#========================================

    def _psi_intercepts(self, ipc):
        jh = len(self.d['psi'][0,:])//2

        return np.abs(self.d['psi'][ipc,:jh] - self.psi0).argmin(),\
               np.abs(self.d['psi'][ipc,jh:] - self.psi0).argmin() + jh

#========================================

    def _load_movie(self, time, slc):
        mvars = 'all'
        d = self._M.get_fields('all', time, slc=slc)

        # This is dumb and should be done in moive
        #if slc is None:
        #    ip,jp = 2*[np.s_[:]]
        #else:
        #    ip,jp = slc[:2]

        #d['xx'] = d['xx'][ip]
        #d['yy'] = d['yy'][jp]

        for q,s in zip([1.,-1.],'ie'):
            for k in 'xx xy xz yy yz zz'.split():
                d['t'+s+k] = d['p'+s+k]/d['n'+s]

            rotate_ten(d,'t'+s, av='')

            for k in 'xyz':
                d['v'+s+k] = q*d['j'+s+k]/d['n'+s]
            
            mag = lambda _x,_y: _x**2 + _y**2

            #d['|b|'] = np.sqrt(reduce(mag, (d['b'+k] for k in'xyz')))
            d['|b|'] = np.sqrt(d['bx']**2 + d['by']**2 + d['bz']**2)

            d['psi'] = gf(calc_psi(d), sigma=self.sig)

        return d

#========================================

    def clean_fig(self):
        import matplotlib.pyplot as plt
        if not hasattr(self, 'fig'):
            self.fig = plt.figure(1,dpi=1200)

        self.fig.clf()

        self.ax = [self.fig.add_subplot(6,2,1+c) for c in range(12)]
        self.fig.set_size_inches(.75*8.5, .75*11.)
        self.fig.subplots_adjust(left=.1,right=.9,bottom=.1,top=.9)

        self.page_counter += 1

#========================================

    def savefig(self, ext='png', dpi=1200):
        fname = '{:03d}_patplots_save.{}'.format(self.page_counter, ext)
        print 'Saving {}...'.format(fname)
        self.fig.tight_layout()
        self.fig.savefig(fname)
        self._created_files += [fname]

#========================================

    def _make_pdf(self, overwrite=0):
        from subprocess import call
        pdfname = date_file_prefix() + 'patplot.pdf'
        print 'Making the pdf...'
        call(['convert']+self._created_files+[pdfname])

        # Lets clean some stuff up for reuse sake
        print 'Removing files...'
        call(['rm','-f'] + self._created_files)

        self._created_files = []
        delattr(self, 'fig')

#========================================

def pat_plots(slc=None, **mvargs):
    """ Make a bunch of figures in the style of Kittypat's work
        
        :param slc: A subslice if you do not want to plot everything
        :type slc: slice or numpy.s_[x, y, z, time]

        :param mvargs: keywork movie arguments to be passed to py3d.Movie
        :type mvargs: kwargs
    """

    pass
