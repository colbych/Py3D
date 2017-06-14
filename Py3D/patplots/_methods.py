import Py3D.sub as ps
import pdb

ctrs = []

#========================================

def plot2D(ax, d, vrs, labs, **kwargs):
    pcm = ps.ims(d, vrs, ax, no_draw=1, **kwargs)
    #pcm = ax.pcolormesh(d['xx'], d['yy'], vrs, **kwargs)

    for c in ctrs:
        ax.plot(*c, color='k',linewidth=.5)

    ax.set_title(labs+':{:.3f} {:.3f}'.format(vrs.min(), vrs.max()), fontsize=6)
    if ax not in ax.get_figure().axes[-2:]:
        ax.set_xlabel('')

    return pcm

#========================================

def gen_conts(d, **ctargs): # d should containe psi
    print 'Generating Contours...'
    import matplotlib.pyplot as plt
    plt.ioff()
    _f,a = plt.subplots(1,1)
    ctargs['colors'] = ctargs.get('colors','k')
    ctargs['linestyles'] = ctargs.get('linestyles', 'solid')
    cts = a.contour(d['xx'], d['yy'], d['psi'].T, **ctargs)
    for c in cts.collections:
        for k in c.get_paths():
            ctrs.append(k.vertices.T)
    
    del(_f)

#========================================

def plot1D(ax, xy, vrs, labs, loc,
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

    if ptargs is not None:
        if type(ptargs) is dict:
            ptargs =3*(ptargs,)
    else:
        ptargs = 3*({},)

    lines = []
    for v,l,pta in zip(vrs, labs, ptargs):
        lines += ax.plot(xy, v, label=l, **pta)

    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)

    ax.set_title('cut @ x = %1.2f'%loc,size=6,loc='right')

    if ldgargs:
        if type(ldgargs) is not dict: ldgargs = {}
        ax.legend(**ldgargs)

#========================================

def psi_intercepts(d, psi0, ipc):
    import numpy as np
    jh = len(d['psi'][0,:])//2
    return np.abs(d['psi'][ipc,:jh] - psi0).argmin(),\
           np.abs(d['psi'][ipc,jh:] - psi0).argmin() + jh

#========================================

def plot_psi_intercepts(ax, d, psi0, ipc):
    yl = ax.get_ylim()
    for kp in psi_intercepts(d, psi0, ipc):
        ax.plot(2*[d['yy'][kp]], yl, 'k--', linewidth=.5)
    ax.set_ylim(yl)

#========================================

def calc_extra_vars(d):
    import numpy as np

    for q,s in zip([1.,-1.],'ie'):
        for k in 'xx xy xz yy yz zz'.split():
            d['t'+s+k] = d['p'+s+k]/d['n'+s]

        ps.rotate_ten(d,'t'+s, av='')

        for k in 'xyz':
            d['v'+s+k] = q*d['j'+s+k]/d['n'+s]
        
        mag = lambda _x,_y: _x**2 + _y**2

        d['|b|'] = np.sqrt(reduce(mag, (d['b'+k] for k in'xyz')))

        d['psi'] = ps.calc_psi(d)

#========================================

def clean_fig():
    import matplotlib.pyplot as plt
    fig = plt.figure(1)
    fig.clf()
    ax = [fig.add_subplot(6,2,1+c) for c in range(12)]
    fig.set_size_inches(.75*8.5, .75*11.)
    return fig,ax

#========================================

def savefig(fig, pgn, ext='png', dpi=300):
    fname = '{:03d}_patplots_save.{}'.format(pgn,ext)
    print 'Saving {}...'.format(fname)
    fig.tight_layout()
    fig.savefig(fname)

