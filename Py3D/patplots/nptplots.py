from scipy.ndimage import gaussian_filter as gf
import _methods

tp = np.s_[500]
ip = np.s_[:8192:2]
jp = np.s_[4096::2]
kp = np.s_[0]

if raw_input('Load data?\n> ') == 'y':
    mvars = 'all'
    M = Movie(0, 'param_asym046')
    d = M.get_fields(mvars, 0, slice=np.s_[tp, kp, jp, ip])
    #d = M.get_fields(mvars, 0)

    dx = d['xx'][1] - d['xx'][0]
    #xy = dict(xx=d['xx'][ip], yy=d['yy'][jp])
    d['xx'] = d['xx'][ip]
    d['yy'] = d['yy'][jp]
    _methods.calc_extra_vars(d)

sig = 3
xy_lim = [[4.5, 6.], [7.3, 8.]]
psi0 = 2.3915547; dpsi = .05; npsi = 10
ctargs = dict(linewidths=.5, 
              levels=arange(psi0 - dpsi*npsi,
                            psi0 + dpsi*npsi,
                            dpsi))

_methods.gen_conts(d, **ctargs)

ldgargs = dict(loc=3,
               ncol=10,
               borderaxespad=0.,
               frameon=0,
               prop={'size':6},
               bbox_to_anchor=(0., 1.02, 1., .102))

#xpcs = [5.025]
xpcs = arange(4.6, 6., .1)
ipcs = [abs(d['xx'] - k).argmin() for k in xpcs]


pgvars = ['ni rho bx ex by ey bz ez |b| jx jy jz'.split(),
          'jix jex jiy jey jiz jez vix vex viy vey viz vez'.split()]

pgvars1D = ['bx by bz |b|', 'ni ne rho',
            'ni ne', 'piyy peyy',
            'ex ey ez', 'jx jy jz',
            'vix viy viz', 'vex vey vez',
            'tixx tiyy tizz', 'texx teyy tezz',
            'tipar tiperp1', 'tepar teperp1']

pgvars1D = [k.split() for k in pgvars1D]

for pgn,pgv in enumerate(pgvars):
    fig,ax = _methods.clean_fig()
    for c,(a,v) in enumerate(zip(ax,pgv)):
        print 'Plotting {}...'.format(v)
        _methods.plot2D(ax=a, 
                        d=d, 
                        vrs=gf(d[v], sigma=sig), 
                        labs=v,
                        cmap='bwr')

        for ipc in ipcs:
            a.plot(d['yy']*0.+d['xx'][ipc], d['yy'], 'k--', linewidth=.5)

        if 'xy_lim' in locals():
            a.set_xlim(xy_lim[0])
            a.set_ylim(xy_lim[1])

    _methods.savefig(fig, pgn)


for pgnc,ipc in enumerate(ipcs):
    fig,ax = _methods.clean_fig()
    for c,(a,vs) in enumerate(zip(ax,pgvars1D)):
        _methods.plot1D(a, d['yy'], 
                        (gf(d[v][ipc], sigma=sig) for v in vs), 
                        labs=vs,
                        loc=d['xx'][ipc],
                        xlim=xy_lim[1],
                        ldgargs=ldgargs)

        _methods.plot_psi_intercepts(a, d, psi0, ipc)
        
    _methods.savefig(fig, pgnc+pgn+1)

