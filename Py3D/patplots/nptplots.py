from scipy.ndimage import gaussian_filter as gf
import numpy as np
import _methods

class PatPlotter(object):
    def __init__(self, **mvargs):
        self._M = Movie(**mvargs)
        self.ctrs = []
        self.page_counter = 0


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


    def make_plots(self, time=None, slc=None)
        if slc is not None:
            slc = slc[::-1] 

        self.d  = self._load_movie()

        # These will need to be arguments
        self.sig = 3
        self.xy_lims = [[4.5, 6.], [7.3, 8.]]

        # This is gonna need a whole thing
        psi0 = 2.3915547; dpsi = .05; npsi = 10

        ctargs = dict(linewidths=.5, 
                      levels=arange(psi0 - dpsi*npsi,
                                    psi0 + dpsi*npsi,
                                    dpsi))

        self._ctrs = self._gen_conts(d, **ctargs)

        ldgargs = dict(loc=3, ncol=10, borderaxespad=0., frameon=0, 
            prop={'size':6}, bbox_to_anchor=(0., 1.02, 1., .102))
        
        # This needs to be automated
        _x,_i = self._gen_xcuts()
        self.xpcs = _x
        self.ipcs = _i

        for page_vars in self.page_vars_2D:
            self.clean_fig()
            self._plot2D(page_vars, **kwargs)
            self._plot2D_cutlines()
            self._plot2D_set_lims()
            self.savefig()

    for pgnc,ipc in enumerate(ipcs):
        self.clean_fig()
        self._plot1D()
        self.plot_psi_intercepts(a, d, psi0, ipc)
            
        self.savefig(fig, pgnc+pgn+1)


#========================================

#========================================

    def _plot2D(self, page_vars, var_labels=None, **kwargs):

        pcm = []

        if var_labels is None: var_labels = page_vars

        for a,v,l in zip(self.ax, page_vars, var_labels):

            print 'Plotting {}...'.format(v)
            vrs = gf(self.d[v], sigma=self.sig)
            pcm += [ps.ims(self.d, vrs, ax=a, no_draw=1, **kwargs)]

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
            for ipc in self._ipcs:
                lines += a.plot(d['yy']*0.+d['xx'][ipc], 
                                d['yy'],
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

    def _gen_xcuts(self):
        xpcs = arange(4.6, 6., .1)
        ipcs = [abs(self.d['xx'] - k).argmin() for k in xpcs]
        return xpcs, ipcs

#========================================

    def _gen_conts(self, **ctargs): # d should containe psi
        d = self.d
        print 'Generating Contours...'
        import matplotlib.pyplot as plt
        plt.ioff()
        _f,a = plt.subplots(1,1)

        ctargs['colors'] = ctargs.get('colors','k')
        ctargs['linestyles'] = ctargs.get('linestyles', 'solid')

        cts = a.contour(d['xx'], d['yy'], d['psi'].T, **ctargs)

        ctrs = []
        for c in cts.collections:
            for k in c.get_paths():
                ctrs.append(k.vertices.T)
        
        del(_f)

        return ctrs

#========================================

    def _plot1D(self, ip, page_vars_1D, var_labels, loc,
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
        """

        if var_labels is None: var_labels = page_vars

        # We will come back it this if we have to
        #if not ptargs: 
        #    ptargs = 3*({},)
        #elif type(ptargs) is dict:
        #    ptargs =3*(ptargs,)

        _xy = self.d[2*cut_dir]
        _lim = _xy[[0,-1]]
        _cut_index = np.s_[ip,:]
        if _cut_dir is 'x': _cut_index = _cut_index[::-1]

        lines = []
        for a,vrs in zip(self.ax, page_vars_1D):

            _sl = []
            for v,l,pwargs in zip(vrs, var_labels, plott_kwargs)]
                _sl += a.plot(_xy, gf(v, sigma=self.sig), label=l)
        
            lines += [_sl]

            a.set_xlim(_lim)
            a.set_title('cut @ x = %1.2f'%loc,size=6,loc='right')

            if ldgargs:
                if type(ldgargs) is not dict: ldgargs = {}
                ax.legend(**ldgargs)

        return lines

#========================================

    def psi_intercepts(self, d, psi0, ipc):
        import numpy as np
        jh = len(d['psi'][0,:])//2
        return np.abs(d['psi'][ipc,:jh] - psi0).argmin(),\
               np.abs(d['psi'][ipc,jh:] - psi0).argmin() + jh

#========================================

    def plot_psi_intercepts(self, ax, d, psi0, ipc):
        yl = ax.get_ylim()
        for kp in psi_intercepts(d, psi0, ipc):
            ax.plot(2*[d['yy'][kp]], yl, 'k--', linewidth=.5)
        ax.set_ylim(yl)

#========================================

    def _load_movie(self, time, slc)
        mvars = 'all'
        d = self._M.get_fields('all', time, slice=slc)

        d['xx'] = d['xx'][ip]
        d['yy'] = d['yy'][jp]

        for q,s in zip([1.,-1.],'ie'):
            for k in 'xx xy xz yy yz zz'.split():
                d['t'+s+k] = d['p'+s+k]/d['n'+s]

            ps.rotate_ten(d,'t'+s, av='')

            for k in 'xyz':
                d['v'+s+k] = q*d['j'+s+k]/d['n'+s]
            
            mag = lambda _x,_y: _x**2 + _y**2

            d['|b|'] = np.sqrt(reduce(mag, (d['b'+k] for k in'xyz')))

            d['psi'] = ps.calc_psi(d)

        return d

#========================================

    def clean_fig(self):
        import matplotlib.pyplot as plt
        try:
            self.fig.clf()
        except AttributeError:
            self.fig = plt.figure(1)

        self.ax = [self.fig.add_subplot(6,2,1+c) for c in range(12)]
        self.fig.set_size_inches(.75*8.5, .75*11.)

        self.page_counter += 1

#========================================

    def savefig(self, pgn, ext='png', dpi=300):
        fname = '{:03d}_patplots_save.{}'.format(page_counter, ext)
        print 'Saving {}...'.format(fname)
        self.fig.tight_layout()
        self.fig.savefig(fname)


def pat_plots( slc=None, **mvargs):
    """ Make a bunch of figures in the style of Kittypat's work
        
        :param slc: A subslice if you do not want to plot everything
        :type slc: slice or numpy.s_[x, y, z, time]

        :param mvargs: keywork movie arguments to be passed to Py3D.Movie
        :type mvargs: kwargs
    """

    pass
