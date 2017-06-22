import numpy as np
import matplotlib.pyplot as plt
from scipy.io.idl import readsav 
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import AutoMinorLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
from Py3D.movie import Movie
from Py3D.dumpID import DumpID

__all__ = ['set_local', 'ims', 'find_xpt', 'var_at', 'ims_subplot',
           'load_movie', 'check_energy_conservation', 'multi_color',
           'show_energy', 'calc_psi', 'readsave', 'date_file_prefix',
           'rs3d', 'rotate_ten']

#======================================================
def set_local(d, loc, overwrite=False):
    """ Sets the contents of a dictory d to the local namespace loc

    Args:
        d (dict): Typicly the field values for a simulation
        loc (dict, namespace): Dictionary that keeps track of the namespace
            you want to load to
        overwrite (bool): If true will overwrite values in loc

    """
    for k in d: #For Key in Simulation Dictionary
        if k not in loc or overwrite:
            loc[k] = d[k]

#======================================================

def ims3D(d,
          k,
          ax=None,
          slice=(2,0),
          extent=None,
          cbar=None,
          cont=None,
          no_draw=None,
          ctargs={},
          **kwargs):
    """ 3D version of the automated imshow (ims)

    Args:
        d (dict): The field values for a simulation, must contain xx, yy & zz
        k (str (key of d) or numpy.array): If k is a str then it will set
            `k = d[k], else k must be a numpy array with the dimensions
            (len(d['xx']), len(d['yy']), len(d['zz'])
        slice (tuple(int, int)): First int is the axis perpendicular to the
            plane of the cut where 012 -> xyz. The second int is the offset
            along the axis. e.g. (2, 20) -> a cut in the X,Y plane taken at
            the 20th index in the Z direction

    The rest are kwargs for ims, see sub.ims

    """

    def try_slice(pv,sl):
        try:
            pv = pv[sl].T
        except IndexError:
            pv = pv.T
        return pv

    old_ax = plt.gca() # Get Current Axis
    if ax is None: 
        ax = old_ax
    else:
        plt.sca(ax)    # Set Current Axis

    if type(k) is str: plt_val = d[k]
    else               : plt_val = k

    if slice[0] == 0:
        plt_val = try_slice(plt_val,np.s_[:,:,slice[1]])
        xlab = r'$Y (d_i)$'; ylab = r'$Z (d_i)$'
        ext = [d['yy'][0], d['yy'][-1], d['zz'][0], d['zz'][-1]]

    elif slice[0] == 1: # We break standered convention because it looks beter
        plt_val = try_slice(plt_val, np.s_[:,slice[1],:])
        plt_val = plt_val.T
        xlab = r'$X (d_i)$'; ylab = r'$Z (d_i)$'
        ext = [d['xx'][0], d['xx'][-1], d['zz'][0], d['zz'][-1]]

    elif slice[0] == 2:
        plt_val = try_slice(plt_val, np.s_[:,:,slice[1]])
        xlab = r'$X (d_i)$'; ylab = r'$Y (d_i)$'
        ext = [d['xx'][0], d['xx'][-1], d['yy'][0], d['yy'][-1]]
    else:
        err_msg = 'slice[0] {} not understood! must be value between 0-2\n'\
                  'where 0 -> Y,Z plane\n'\
                  '      1 -> Z,x plane\n'\
                  'and   2 -> X,Y plane'
        print err_msg.format(slice[0])
        raise IOError()
    
    return_tuple = _ims(d,plt_val,xlab,ylab,ext,ax,extent,cbar,
                         cont,no_draw,ctargs,**kwargs)
    plt.sca(old_ax)
    return return_tuple 

#======================================================

def ims(d,
        k,
        ax=None,
        extent=None,
        cbar=None,
        cont=None,
        no_draw=None,
        ctargs={},
        **kwargs):
    """
    A wrapper function for imshow to do most tedious stuff for P3D simulations
    
    Args:
        d (dict): A dictionary with relevent simualtion information.
            d must contain xx and yy so it will know the dimensions to plot
        k (str or np.array): Either a str of a varible contained within d or 
        a 2D numpy array of size (len(d['xx']), len(d['yy'])) that will be 
        plotted.
        cbar (bool): If true, then auto generate a colorbar
        cont (bool): If true, then auto generate contours
        no_draw (bool): If ture, do not call matplotlib.pylab.draw(), can
            speed up the plotting process
        ctargs (dict): A dictonary to pass extra argumens to the contour
            call, so you can add more lines or set the elvels explicitly.

    """

    old_ax = plt.gca() # Get Current Axis
    if ax is None: 
        ax = old_ax
    else:
        plt.sca(ax)    # Set Current Axis

    if type(k) is str: plt_val = d[k]
    else               : plt_val = k
    plt_val = plt_val.T


    xlab = r'$X (d_i)$'
    ylab = r'$Y (d_i)$'

# Use the dict values of xx and yy to set extent
    ext = [d['xx'][0],
           d['xx'][-1],
           d['yy'][0],
           d['yy'][-1]]

    return_tuple = _ims(d,plt_val,xlab,ylab,ext,ax,extent,cbar,
                         cont,no_draw,ctargs,**kwargs)

    plt.sca(old_ax)

    return return_tuple 

#======================================================

def _ims(d,
         plt_val,
         xlab,
         ylab,
         ext,
         ax,
         extent,
         cbar,
         cont,
         no_draw,
         ctargs,
         **kwargs):


    if kwargs.has_key('cmap'): cmap=kwargs.pop('cmap')
    else:                      cmap='PuOr'

    im = ax.imshow(plt_val,
                   origin='low',
                   extent=ext,
                   cmap=cmap,            # I just love this color map
                   **kwargs)

    if extent is not None:
        ax.set_xlim(extent[:2])
        ax.set_ylim(extent[2:])

    ax.autoscale(False)
    ax.set_xlabel(xlab,size=8)
    ax.set_ylabel(ylab,size=8)

    ax.xaxis.set_tick_params(which='both',labelsize=8)
    ax.yaxis.set_tick_params(which='both',labelsize=8)
    plt.minorticks_on()

    return_tup = im,
    # Code to implement for a cbar
    if cbar:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", "3%", pad="1.5%")
        plt.colorbar(im, cax=cax)

        cax.xaxis.set_tick_params(which='both',labelsize=8)
        cax.yaxis.set_tick_params(which='both',labelsize=8)

        return_tup += (cax,)

    if cont or ctargs:
        if 'psi' in d:
            psi = d['psi']
        else:
            psi = calc_psi(d)
        if 'colors' not in ctargs: ctargs['colors'] = 'k'
        if 'linestyles' not in ctargs: ctargs['linestyles'] = 'solid'
       
        cts = ax.contour(d['xx'],d['yy'],psi.T,**ctargs)

        return_tup += (cts,)

    if not no_draw:
        plt.draw()

    if len(return_tup) == 1:
        return return_tup[0]
    else:
        return return_tup

#======================================================

def find_xpt(d):
    psi = calc_psi(d)
    jp = int(np.round(len(d['yy'])/2.))
    yp = d['yy'][jp]

    if 'bxav' in d: av = 'av'
    else: av = ''
    lBm = np.mean(d['bx'+av][:jp,:])
    uBm = np.mean(d['bx'+av][jp:,:])

    if lBm > uBm: #Upper
        ip = psi[jp,:].argmax()
        print 'Finding upper max'
    else:         #Lower
        ip = psi[jp,:].argmin()
        print 'Finding lower min'
    xp = d['xx'][ip]

    return ip,jp,xp,yp

#======================================================
    
def ims_subplot(d, var, ax, window, **kwargs):
    """ ims, but for a sub region, plots faster and get auto max and min
        colors right
    
    Args:
        d (dict): Typicly the field values for a simulation
        k (str (key of d) or numpy.array): the string of field that you want
            to plot.
        ax (matplotlib.pyplot.axis): MPL axis to plot to
        window (4*[float]): List of corners of the window to plot
    """

    ix = lambda x,xx: np.abs(xx - x).argmin()
    iwd = [ix(w,d[v]) for w,v in zip(window,['xx','yy','xx','yy'])]
    
    av = 'av' if 'bxav' in d else ''

    td = {'xx':d['xx'][iwd[0]:iwd[2]],
          'yy':d['yy'][iwd[1]:iwd[3]],
          'bx'+av:d['bx'+av][iwd[1]:iwd[3], iwd[0]:iwd[2]],
          'by'+av:d['by'+av][iwd[1]:iwd[3], iwd[0]:iwd[2]],
          'bz'+av:d['bz'+av][iwd[1]:iwd[3], iwd[0]:iwd[2]]}

    rst = ims(td, var[iwd[1]:iwd[3], iwd[0]:iwd[2]], ax, **kwargs)

    return rst 

#======================================================

def shift_to_xpt_frame(d):
    ip,jp,xp,yp = find_xpt(d)
    d['xx'] = d['xx'] - xp  
    d['yy'] = d['yy'] - yp  

#======================================================

def var_at(fdic, key, r0, ordflg='idl'):
    delx = fdic['xx'][1] - fdic['xx'][0] 
    dely = fdic['yy'][1] - fdic['yy'][0] 

    xind = (np.floor((r0[0]-delx/2.0)/delx)).astype(int)
    yind = (np.floor((r0[1] - fdic['yy'][0])/dely)).astype(int)

    wx = (r0[0]-delx/2.0)%delx
    wy = (r0[1]-fdic['yy'][0])%dely

    if ordflg =='idl':
        var = wx     *wy     *fdic[key][yind    ,xind    ] + \
              (1.-wx)*wy     *fdic[key][yind    ,(xind+1)] + \
              wx     *(1.-wy)*fdic[key][(yind+1),xind    ] + \
              (1.-wx)*(1.-wy)*fdic[key][(yind+1),(xind+1)] 
    else:
        var = wx     *wy     *fdic[key][xind    ,yind] + \
              (1.-wx)*wy     *fdic[key][(xind+1),yind] + \
              wx     *(1.-wy)*fdic[key][xind    ,(yind+1)] + \
              (1.-wx)*(1.-wy)*fdic[key][(xind+1),(yind+1)] 

    return var

#======================================================

def load_movie(num=None,
               param=None,
               path='./',
               vrs='all',
               time=None,
               slc=None,
               name_style='p3d'):

    """ Parameters
        ----------
        num : int
            Moving number to load. If None it will ask
        param : str
            name of param file, If None it will ask.
        path : str
            where movie files are. (Assumese in local dir)
        vrs (vars) : str or array of strs
            what varibles to load ['bx', 'by', 'bz', ..],
            Assumes that you want to load 'all'
        slc (slice) : slice
            A slice from the movie. It is easiest to just pass 
            np.s_[x, y, z, t], where x,y,z,t are the the range
            over the coresponding axes
        time : int
            what time to load the move from. If None it will ask
        name_style : str ('p3d' or 'tulasi' or 'unfinished')
            Tulasi's version of the code has different nameing convention
            for movie files, so you have to used the UnfinsihedMovie
            object to get this to work.
    """

    if name_style is 'p3d':
        return Movie(num,param,path).get_fields(vrs,time,slc)
    elif name_style.lower() in ['tulasi', 'unfinished'] :
        from Py3D.movie import UnfinishedMovie
        return UnfinishedMovie(param,path).get_fields(vrs,time,slc)



#======================================================

def gen_distro(species,
               r=[1.,1.],
               dx=[.5,.5],
               par=False,
               **vdargs):
    raise NotImplementedError()
    
#======================================================

def check_energy_conservation(mov_num=0,
                              init_time=0, 
                              final_time=-1,
                              ims_var='pexx',
                              use_UFM=False):

    """ A simpile method that plots how well the total energy is conserved
        Parameters
        ----------
        mov_num 

    """
    from scipy.ndimage import gaussian_filter as gf
    import glob

    vs = 'pexx peyy pezz pixx piyy pizz ne ni'.split()
    if ims_var not in vs: vs.append(ims_var)

    print 'Loading Energy from p3d.stdout.000...'
    try:
        engs = show_energy('p3d.stdout.000').astype(float)
    except:
        print 'p3d.stdout.000 file not found!'
        return None

    try:
        param_file = glob.glob('param*')[0]

        if use_UFM:
            from Py3D.movie import UnfinishedMovie
            print 'Loading temperatures from movies...'
            M = UnfinishedMovie(param_file)
        else:
            print 'Loading temperatures from movie.???.000...'
            M = Movie(0,param_file)
        slc = (2,0) if M.param['nz']*M.param['pez'] > 1 else None

        if final_time < 0 : final_time += M.ntimes #allows for negetive indexing

        di = M.get_fields(vs, init_time, slice=slc)
        df = M.get_fields(vs, final_time, slice=slc)

    except:
        print 'Moive did not load properly! Exiting!!!'
        return None
    
    fig = plt.figure(1)
    fig.clf()
    fig.set_size_inches(12.8, 8.2)
    ax = [fig.add_subplot(331+c) for c in range(9)]
    
    tvc = np.arange(engs.shape[0])*M.param['dt']

    ### subplot 1,1
    a =ax[0]
    for eng,l in zip(engs.T, 'B KE Tot'.split()):
        a.plot(tvc,eng,label=l)
    a.set_xlabel('time'); a.set_ylabel('Energy')

    ### subplot 1,2
    a = ax[1]
    for eng,l in zip(engs.T, 'B KE Tot'.split()):
        a.plot(tvc,100.*(eng - eng[0])/eng[0],label=l)
    a.set_xlabel('time'); a.set_ylabel('% Energy change')

    a = ax[2]
    ims(df, ims_var, ax=a, no_draw=1)

    jp = len(di['yy'])/2
    sig = 3.

    ### subplot 2,1
    a = ax[3]
    for k in 'xyz':
        lab = 'Te'+2*k+', t='+str(tvc[0])
        a.plot(di['xx'], (di['pe'+2*k]/di['ne'])[:,0], label=lab)
        lab = 'Te'+2*k+', t='+str(tvc[-1])
        a.plot(df['xx'], (df['pe'+2*k]/df['ne'])[:,0], label=lab)
    a.set_xlabel('X di'); a.set_ylabel('Te')
    a.legend(loc='best',ncol=3,fontsize=6)

    ### subplot 2,2
    a = ax[4]
    for k in 'xyz':
        lab = 'Te'+2*k+', t='+str(tvc[0])
        a.plot(di['xx'], (di['pe'+2*k]/di['ne'])[:,jp], label=lab)
        lab = 'Te'+2*k+', t='+str(tvc[-1])
        a.plot(df['xx'], (df['pe'+2*k]/df['ne'])[:,jp], label=lab)
    a.set_xlabel('X di'); a.set_ylabel('Te')

    ### subplot 2,3
    a = ax[5]
    for k in 'xyz':
        lab = 'Te'+2*k+', t='+str(tvc[0])
        a.plot(di['xx'], gf((di['pe'+2*k]/di['ne'])[:,0], sigma=sig), label=lab)
        lab = 'Te'+2*k+', t='+str(tvc[-1])
        a.plot(df['xx'], gf((df['pe'+2*k]/df['ne'])[:,0], sigma=sig), label=lab)
    a.set_xlabel('X di'); a.set_ylabel('Te')

    ### subplot 3,1
    a = ax[6]
    for k in 'xyz':
        lab = 'Ti'+2*k+', t='+str(tvc[0])
        a.plot(di['xx'], (di['pi'+2*k]/di['ni'])[:,0], label=lab)
        lab = 'Ti'+2*k+', t='+str(tvc[-1])
        a.plot(df['xx'], (df['pi'+2*k]/df['ni'])[:,0], label=lab)
    a.set_xlabel('X di'); a.set_ylabel('Ti')

    ### subplot 3,2
    a = ax[7]
    for k in 'xyz':
        lab = 'Ti'+2*k+', t='+str(tvc[0])
        a.plot(di['xx'], (di['pi'+2*k]/di['ni'])[:,jp], label=lab)
        lab = 'Ti'+2*k+', t='+str(tvc[-1])
        a.plot(df['xx'], (df['pi'+2*k]/df['ni'])[:,jp], label=lab)
    a.set_xlabel('X di'); a.set_ylabel('Ti')

    ### subplot 3,3
    a = ax[8]
    for k in 'xyz':
        lab = 'Ti'+2*k+', t='+str(tvc[0])
        a.plot(di['xx'], gf((di['pi'+2*k]/di['ni'])[:,0], sigma=sig), label=lab)
        lab = 'Ti'+2*k+', t='+str(tvc[-1])
        a.plot(df['xx'], gf((df['pi'+2*k]/df['ni'])[:,0], sigma=sig), label=lab)
    a.set_xlabel('X di'); a.set_ylabel('Ti')

    plt.tight_layout()


#======================================================

def multi_color(slice=None, draw=False):

    """ A method for Mike!. It coppies his multi gray
        IDL code.

        Parameters
        ----------
        slice : None || (axis, offset)
            Only relevent for 3D data. Which slice (plane) to plot.
            axis : int 0,1,2
                Which plane do you want to see? 
                0 -> (y,z), 1 -> (x,z), 2 -< (x,y)
            offset : float
                What value do you want to use in the 3rd direction?

        draw : bool
            if True it will plot and draw every subplot in real time
            Note: if set to True it is VERY slow
    """
    M,t,fig,istate = _movie_start_plot()

    # This is going to sound crazy but there seems to be a BIG time
    # difference in loading the movie files depending on if the slice
    # is a list or a tuple (like a factor of 10 faster for tuple)
    if slice is None: 
        slice=(2,0)
    else:
        xyz = M._get_xyz_vectors()
        k = 2*('xyz'[slice[0]])
        s2 = np.argmin(np.abs(xyz[k] - slice[1]))
        slice = (slice[0], s2)


    print 'Slice = ',slice
    print 'Making subplots...'
    ax = [fig.add_subplot(6,5,c+1) for c in range(6*5)]
    im = []
    for a,k in zip(ax,M.movie_vars):

        print 'loading ',k
        d = M.get_fields(k, time=t, slice=slice)

        print 'plotting ',k
        ttl = k
        if M.param['pez']*M.param['nz'] > 1:
            im.append(ims3D(d,k,a, no_draw=not draw, slice=slice))
            ttl+= ', {}={}'.format('xyz'[slice[0]], slice[1])
            a.set_title(ttl,size=8)

        else:
            im.append(ims(d,k,a, no_draw=not draw))
            a.set_title(ttl,size=8)
    
    _movie_end_plot(istate)
    return ax,im

#======================================================

def three_plane(v, r0=[0., 0., 0.], **imsargs):
    """ A method to look at 3 planes that intersect a point

        Parameters
        ----------
        v : str
            (v)ariable to plot, there must be a coresponding movie file
        r0 : (x, y, z)
            x : float
                x position of cut in the (y,z) plane
            y : float
                y position of cut in the (x,z) plane
            z : float
                z position of cut in the (x,y) plane
        imsargs : key word argumetns for Py3D.sub.ims3D
    """
    
    M,t,fig,istate = _movie_start_plot()
    
    # Note this is a sloopy way of converting r0 -> ind0 this takes O(N) 
    # time where N is the length of a given dimentison. We can easly do 
    # this in O(1) time
    xyz = M._get_xyz_vectors()
    ind0 = [np.argmin(np.abs(xyz[k+k] - r)) for r,k in zip(r0,'xyz')]

    print 'Making subplots...'
    ax = [fig.add_subplot(2,2,c+1) for c in [0,3,2]]
    im = []
   
    for c in range(3):
        print 'loading x,y slice of {} @ kp_z={}'.format(v,ind0[2])

        slice = ((c+2)%3,ind0[(c+2)%3])

        d = M.get_fields(v, time=t, slice=slice)
        im.append(ims3D(d,v,ax[c], no_draw=True, slice=slice, **imsargs))
        ttl = '{}={}'.format('xyz'[slice[0]], r0[(c+2)%3])
        ax[c].set_title(ttl,size=8)

    _x2 = max([xyz[k][-1] for k in xyz])
    _x1 = min([xyz[k][0]  for k in xyz])
    _x = np.linspace(_x1,_x2,200)
    for a,y,x in zip(ax, [0,1,0], [1,2,2]):
        a.plot(0.*_x + r0[0], _x, 'k--')
        a.plot(_x, 0.*_x + r0[1], 'k--')
    
    if 'vmin' not in imsargs and 'vmax' not in imsargs:
        def gcm(_im):
            if type(_im) is tuple:
                _im = _im[0]
            return _im

        cmax = max([gcm(_im).get_clim()[1] for _im in im])
        cmin = min([gcm(_im).get_clim()[0] for _im in im])

        for _im in im:
            gcm(_im).set_clim(cmin, cmax)

    _movie_end_plot(istate)

    return ax,im
        

#======================================================

def _movie_start_plot():
    M = Movie()
    t = raw_input('Enter time between {}-{}: '.format(0,M.ntimes-1))
    t = int(t)

    print 'Getting Fgiure...'
    fig = plt.gcf()
    fig.clf()
    istate = plt.isinteractive
    plt.ioff()

    return M,t,fig,istate

def _movie_end_plot(istate):
    if istate:
        plt.ion()
    plt.tight_layout()
    plt.draw()

#======================================================


def show_energy(fname=None):
    """ Grabs the energy values from a p3d.stdout file and returns them
        as a numpy array.
    
    Args:
        fname (str, optional): Name of the p3d.stdout file to grab.
            If None it will ask.
    """
    if fname is None:
        fname = raw_input('Enter p3d.stdout file: ')

    f = open(fname, 'r')
    eng = []
    for lines in f:
        if lines.find('ENERGY') > -1 and lines.find('ENERGY:') < 0:
            eng.append(lines.split()[1:4])
    f.close()

    return np.array(eng).astype(float)

#======================================================

def calc_psi(d):
    """ Calculated the magnetic scaler potential for a 2D simulation

    Args:
        d (dict): Dictionary containing the fields of the simulation
            d must contain bx, by, xx and yy

    Retruns:
        psi (numpy.array(len(d['xx'], len(d['yy']))) ): Magnetic scaler
            potential
    """

    if 'bxav' in d and 'byav' in d:
        bx = d['bxav']
        by = d['byav']
    else:
        bx = d['bx']
        by = d['by']

    psi = 0.0*bx
    psi[0,1:] = np.cumsum(bx[0,1:])*(d['yy'][2] - d['yy'][1])
    psi[1:,:] = psi[0,:] - np.cumsum(by[1:,:], axis=0)*(d['xx'][2] - d['xx'][1])
    return psi

#======================================================

def plot_line(itcpt,dir='y',ax=None,**kwargs):
    if ax is None:
        ax = plt.gca()
    
    xarr = array(ax.get_xlim())
    yarr = array(ax.get_ylim())

    if dir == 'x':
        yarr = yarr*0.0 + itcpt
    elif dir == 'y':
        xarr = xarr*0.0 + itcpt
    else:
        print 'I dont understand what direction ' + str(dir) + \
              'is! Nothing plotted.'
        return None

    ax.plot(xarr,yarr,**kwargs)

#======================================================
# Method hold over from p3d_run, needs to be revised!
#
#def avg_movie(fname=None, 
#              param=None,
#              mov=None,
#              ntimes=None,
#              save_full=None):
#    
#    way = 'slow'
#
#    if fname is None:
#        fname = raw_input('Enter Save file name base: ')
#
#    if way == 'fast':
#        CC = p3d_run('local',param=param)
#        print 'Loading time %i'%0
#        CR = CC.load_movie('all',0)
#
#        ntimes = CC.movie.num_of_times
#
#        keys = CR.keys()
#        keys.pop(keys.index('xx'))
#        keys.pop(keys.index('yy'))
#
#        t = time.time()
#        for k in keys:
#            _ = CC.load_movie(k,range(1,ntimes),mov)
#            CR[k] += np.sum(_[k],axis=0) # sum along time
#            CR[k+'av'] = CR[k]/1.0/ntimes
#            CR[k] = _[k][-1,:,:]
#
#        print 'TOTAL TIME: %f'%(time.time() - t)
#
#    else:
#    ## First way I tried, may be slow?
#        CC = p3d_run('local',param=param)
#        print 'Loading time %i'%0
#        CR = CC.load_movie('all',0,mov)
#        
#        if ntimes is None:
#            ntimes = CC.movie.num_of_times
#
#        keys = CR.keys()
#        keys.pop(keys.index('xx'))
#        keys.pop(keys.index('yy'))
#
#        t = time.time()
#        for cosa in range(1,ntimes):
#            print '\n==================\n' \
#                    'Loading time %i' \
#                  '\n==================\n'%cosa
#            _ = CC.load_movie('all',cosa)
#            for k in keys:
#                CR[k] += _[k]
#                if cosa == ntimes -1:
#                    CR[k+'av'] = 1.0*CR[k]/ntimes
#                    CR[k] = _[k]
#
#        rotate_ten(CR,'pi')
#        rotate_ten(CR,'pe')
#        CR['tiparav'] = CR['piparav']/CR['niav']
#        CR['tiperp1av'] = CR['piperp1av']/CR['niav']
#        CR['tiperp2av'] = CR['piperp2av']/CR['niav']
#        CR['teparav'] = CR['peparav']/CR['neav']
#        CR['teperp1av'] = CR['peperp1av']/CR['neav']
#        CR['teperp2av'] = CR['peperp2av']/CR['neav']
#
#        print 'TOTAL TIME: %f'%(time.time() - t)
#
#        CRL = {}
#        CRU = {}
#        ylen = len(CR['yy'])
#        for k in CR:
#            if k != 'yy' and k != 'xx':
#                CRL[k] = np.squeeze(CR[k])[:ylen/2,:]
#                CRU[k] = np.squeeze(CR[k])[ylen/2:,:]
#            elif k == 'yy':
#                CRL[k] = CR[k][:ylen/2]
#                CRU[k] = CR[k][ylen/2:]
#            else:
#                CRL[k] = CR[k]
#                CRU[k] = CR[k]
#        
#        print 'Saving lower data...'
#        np.save(fname+'_lower',CRL)
#        print 'Saving upper data...'
#        np.save(fname+'_upper',CRU)
#
#        if save_full:
#            np.save(fname,CR)
#
#    return CR


#======================================================

def roll_run(d, sx=None):
    """ Roll every variable in a simulation (d)
        in the x direction by length in indexspace (sx)

        Warning: Old method from p3dthon use with catution
    """
    klst = ['rho','jx','jy','jz','bx','by','bz',
            'ex','ey','ez','ne','jex','jey','jez',
            'pexx','peyy','pezz','pexy','peyz','pexz',
            'ni','jix','jiy','jiz','vix','viy','viz',
            'vex','vey','vez','dene','deni',
            'pixx','piyy','pizz','pixy','piyz','pixz',
            'tepar','teperp1','tipar','tiperp1']

    kavlst = [k+'av' for k in klst]
    if sx is None:
        if d['yy'][0] < 1.0: 
            print 'Gonna roll RIGHT!!!!!!!!!!!!'
            sx = -1*np.size(d['xx'])/4
        else: 
            print 'Gonna roll LEFT!!!!!!!!!!!!'
            sx = np.size(d['xx'])/4
    
    for key in d.keys():
# Old way not super smart
#        if key.rfind('av') == len(key)-2 and len(key) > 2:
#            print 'Rolling ',key
#            d[key] = np.roll(d[key],sx,axis=1)
# New way a little bit smarter
        if key in klst or key in kavlst:
            print 'Rolling ',key
            d[key] = np.roll(d[key],sx,axis=1)

#======================================================

def readsave(restore_fname):
    """ read an idl .sav(.dat) file or a python npy file """
    if restore_fname[restore_fname.rfind('.'):] == '.npy':
        d = np.load(restore_fname).all()
        for v in ['ni', 'ne', 'niav', 'neav']:
            if v in d:
                d['de'+v] = d[v]
        return d
    else:
        return readsav(restore_fname)

#======================================================

def date_file_prefix():
    """ returns a string with the current date """
    import datetime
    return datetime.date.today().strftime('%Y.%m.%d.')

#======================================================

def rs3d(arr):
    """ Reshape an array as a 3D array (for Tulasi's stupid code)"""
    return arr.reshape(arr.shape + (1,))

#======================================================

def rotate_ten(d,
               var='pi',
               av='av',
               overwrite=False,
               full_rotate=False):

    if var+'par'+av in d and not overwrite:
        print 'Warning: {} was found in the'.format(var+'par'+av) +\
              'restored data: nothing will be rotated!!!!'
        pass

        
    elif full_rotate:
# e1 -> \hat{B} 
# e2 -> \hat{ExB}
# e3 -> \hat{Bx(ExB)}
        e1 = np.array([d['bxav'],
                       d['byav'],
                       d['bzav']])

        e2 = np.cross(np.array([d['exav'],
                                d['eyav'],
                                d['ezav']]), e1,axis=0)

        e1 = e1/np.sqrt(np.sum(e1**2,axis=0))
        e2 = e2/np.sqrt(np.sum(e2**2,axis=0))
        e3 = np.cross(e1,e2,axis=0)

        T = np.array([[d[var+'xx'+av],d[var+'xy'+av],d[var+'xz'+av]],
                      [d[var+'xy'+av],d[var+'yy'+av],d[var+'yz'+av]],
                      [d[var+'xz'+av],d[var+'yz'+av],d[var+'zz'+av]]])

        Te1 = np.array([np.sum(T[0,:,:,:]*e1,axis=0),
                        np.sum(T[1,:,:,:]*e1,axis=0),
                        np.sum(T[2,:,:,:]*e1,axis=0)])
        Te2 = np.array([np.sum(T[0,:,:,:]*e2,axis=0),
                        np.sum(T[1,:,:,:]*e2,axis=0),
                        np.sum(T[2,:,:,:]*e2,axis=0)])
        Te3 = np.array([np.sum(T[0,:,:,:]*e3,axis=0),
                        np.sum(T[1,:,:,:]*e3,axis=0),
                        np.sum(T[2,:,:,:]*e3,axis=0)])

        Tpar  = np.sum(e1*Te1,axis=0)
        Tperp = (T[0,0,:,:] + T[1,1,:,:] + T[2,2,:,:] - Tpar)/2.0

# The diaganal is easy, we pick perp1 = perp2
        d[var+'par'+av]   = Tpar
        d[var+'perp1'+av] = Tperp
        d[var+'perp2'+av] = Tperp

        a = np.sum(e2*Te1,axis=0)
        b = np.sum(e3*Te1,axis=0)
        c = np.sum(e2*Te2,axis=0)
        d = np.sum(e3*Te2,axis=0)
        e = np.sum(e3*Te3,axis=0)
        
        x = (c - e)/d/2.
        ct = np.sqrt(1. + 1./np.sqrt((1.+ (x)**2)))/np.sqrt(2)
        st = x/(np.sqrt(2.)*np.sqrt(x**2.+1.)*\
                np.sqrt(1./np.sqrt(x**2.+1.)+1.))

        T11 = Tpar
        T12 = a*ct - b*st
        T13 = a*st + b*ct
        T22 = c*ct**2 + e*st**2 - 2.*d*st*ct
        T23 = (c - e)*st*ct + d*(ct**2 - st**2)
        T33 = c*st**2 + e*ct**2 + 2.*d*st*ct

        d[var+'11'+av] = T11
        d[var+'12'+av] = T12
        d[var+'13'+av] = T13
        d[var+'22'+av] = T22
        d[var+'23'+av] = T23
        d[var+'33'+av] = T23

        # Now Agyrotropy code!
        d[var+'agy'+av] = 2.*np.sqrt((T12**2 + T13**2 + T23**2)) \
                                 /(T11 + T22 + T33) 
    else:
# This was the old way, and it was very simple

        bmag = np.sqrt( d['bx'+av]**2+
                        d['by'+av]**2+
                        d['bz'+av]**2)
        bbx,bby,bbz = (d[k+av]/bmag for k in 'bx by bz'.split())

        d[var+'par'+av] = (bbx*(bbx*d[var+'xx'+av] + 
                                 bby*d[var+'xy'+av] + 
                                 bbz*d[var+'xz'+av])+
                            bby*(bbx*d[var+'xy'+av] + 
                                 bby*d[var+'yy'+av] + 
                                 bbz*d[var+'yz'+av])+
                            bbz*(bbx*d[var+'xz'+av] + 
                                 bby*d[var+'yz'+av] + 
                                 bbz*d[var+'zz'+av]))

        d[var+'perp1'+av] = (d[var+'xx'+av] + 
                              d[var+'yy'+av] +
                              d[var+'zz'+av] - 
                              d[var+'par'+av])/2.

        d[var+'perp2'+av] = d[var+'perp1'+av]

#======================================================

def calc_pdf(ar, pdf_min=None, pdf_max=None, weight=100, inc=0, ax=0):
    """
    This is Tulasi's code it allegedly works. I have not extensivly
    tested it myself use at your own risk.

    """
    
    if len(ar) == 0:
        print 'No array provided! Exiting!'
        return
    if pdf_min is None:
        pdf_min = ar.min()

    if pdf_max is None:
        pdf_max=ar.max()

    # If PDF of increment, then increment the array
    if inc > 0:
        ar = ar - np.roll(ar,inc,axis=ax)
    
    # Find the total length of data set
    # Also this used to use the operator module, but I didn't like that
    arsize = reduce(lambda x,y: x*y, np.shape(ar),1)

    # Find the RMS value of data set and normalize to it.
    rmsval = np.sqrt(np.mean(ar**2))
    if rmsval != 0:
        ar = ar/rmsval

    # Reshape the array to 1D & sort it.
    arr=np.reshape(ar,arsize)
    np.ndarray.sort(arr,kind='heapsort')

    # Empty arrays for output
    bins=int(arsize/weight); pdf=np.zeros(bins); binvals=np.zeros(bins)

    # Fill the bins 
    for i in range(bins):
        start=i*weight
        binvals[i] = np.mean(arr[start:start+weight])
        pdf[i] = weight/(arr[start:start+weight].max()-arr[start:start+weight].min())
    pdf = pdf/arsize
    return binvals,pdf
