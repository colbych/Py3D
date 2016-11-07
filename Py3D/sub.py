import time
import operator
import pdb
import numpy as np
import matplotlib.pyplot as plt
from scipy.io.idl import readsav 
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import AutoMinorLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
from Py3D.movie import Movie
from Py3D.dumpID import DumpID

#======================================================
def set_local(d,loc,overwrite=False):
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
    A wrapper function for imshow to do most tedious stuff for
    my simulations

    Parameters
    ==========
    d : Simulation Dictionary
        dictionary will relevent simualtion information
        all d really eeds to contain is xx and yy so it
        will know the dimensions

    k : str or 2D Array
        Either a str of a varible contained within d or 
        a 2D array that will be plotted. I usaly usely
        ues the str for 'bz' but the var for vix= jix/ni

    cbar: bool
        1 or 0 on weather to automaticly generage a colorbar

    cont: bool
        1 or 0 on weather to automaticly add field lines to
        the plot

    no_draw : bool
        Set True to skip the renering processees. Makes ims
        a little bit faster

    ctargs : dict
        a dictonary to pass extra argumens to the contour
        call, so you can add more lines or set the elvels
        explicitly.

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

    if cont:
        if 'psi' in d:
            psi = d['psi']
        else:
            psi = calc_psi(d)
        if 'colors' not in ctargs: ctargs['colors'] = 'k'
        if 'linestyles' not in ctargs: ctargs['linestyles'] = 'solid'
        
        cts = ax.contour(d['xx'],d['yy'],psi,**ctargs)

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
    
def ims_subplot(d,var,ax,window,**kwargs):
    ix = lambda x,xx: np.abs(xx - x).argmin()
    iwd = [ix(w,d[v]) for w,v in zip(window,['xx','yy','xx','yy'])]
    
    if 'bxav' in d: av = 'av'
    else: av = ''
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

def var_at(fdic,key,r0,ordflg='idl'):
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
               vars='all',
               time=None):

    """ Parameters
        ----------
        num : int
            Moving number to load. If None it will ask
        param : str
            name of param file, If None it will ask.
        path : str
            where movie files are. (Assumese in local dir)
        vars : str or array of strs
            what varibles to load ['bx', 'by', 'bz', ..],
            Assumes that you want to load 'all'
        time : int
            what time to load the move from. If None it will ask
    """
    return Movie(num,param,path).get_fields(vars,time)

#======================================================

def gen_distro(species,
               r=[1.,1.],
               dx=[.5,.5],
               par=False,
               **vdargs):
    raise NotImplementedError()
    

#======================================================

def multi_color(slice=None, draw=False):

    """ A method for Mike!. It coppies his multi gray
        IDL code.
    """
    plt_plane = 2
    plt_offset = 0

    M = Movie()
    t = raw_input('Enter time between {}-{}: '.format(0,M.ntimes-1))
    t = int(t)
    
    print 'Getting Fgiure...'
    fig = plt.gcf()
    fig.clf()
    plt.ioff()

    print 'Making subplots...'
    ax = [fig.add_subplot(6,5,c+1) for c in range(6*5)]
    for a,k in zip(ax,M.movie_vars):

        print 'loading ',k
        d = M.get_fields(k, time=t, slice=slice)

        print 'plotting ',k
        ttl = k
        if M.param['pez']*M.param['nz'] > 1:
            ims3D(d,k,a, no_draw=not draw, slice=slice)
            ttl+= ', {}={}'.format('xyz'[slice[0]], slice[1])
            a.set_title(ttl,size=8)

        else:
            ims(d,k,a, no_draw=not draw)
            a.set_title(ttl,size=8)
    
    plt.ion()
    plt.draw()
    plt.tight_layout()

        

#======================================================


def show_energy(fname=None):
    if fname is None:
        fname = raw_input('Enter p3d.stdout file: ')

    f = open(fname, 'r')
    eng = []
    for lines in f:
        if lines.find('ENERGY') > -1 and lines.find('ENERGY:') < 0:
            eng.append(lines.split()[1:4])
    f.close()

    return np.array(eng)

#======================================================

def calc_psi(d):
# Calculating Psi                                                                                   
    if 'bxav' in d and 'byav' in d:
        bx = d['bxav']
        by = d['byav']
    else:
        bx = d['bx']
        by = d['by']

    #pdb.set_trace()
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

def avg_movie(fname=None, 
              param=None,
              mov=None,
              ntimes=None,
              save_full=None):
    
    way = 'slow'

    if fname is None:
        fname = raw_input('Enter Save file name base: ')

    if way == 'fast':
        CC = p3d_run('local',param=param)
        print 'Loading time %i'%0
        CR = CC.load_movie('all',0)

        ntimes = CC.movie.num_of_times

        keys = CR.keys()
        keys.pop(keys.index('xx'))
        keys.pop(keys.index('yy'))

        t = time.time()
        for k in keys:
            _ = CC.load_movie(k,range(1,ntimes),mov)
            CR[k] += np.sum(_[k],axis=0) # sum along time
            CR[k+'av'] = CR[k]/1.0/ntimes
            CR[k] = _[k][-1,:,:]

        print 'TOTAL TIME: %f'%(time.time() - t)

    else:
    ## First way I tried, may be slow?
        CC = p3d_run('local',param=param)
        print 'Loading time %i'%0
        CR = CC.load_movie('all',0,mov)
        
        if ntimes is None:
            ntimes = CC.movie.num_of_times

        keys = CR.keys()
        keys.pop(keys.index('xx'))
        keys.pop(keys.index('yy'))

        t = time.time()
        for cosa in range(1,ntimes):
            print '\n==================\n' \
                    'Loading time %i' \
                  '\n==================\n'%cosa
            _ = CC.load_movie('all',cosa)
            for k in keys:
                CR[k] += _[k]
                if cosa == ntimes -1:
                    CR[k+'av'] = 1.0*CR[k]/ntimes
                    CR[k] = _[k]

        rotate_ten(CR,'pi')
        rotate_ten(CR,'pe')
        CR['tiparav'] = CR['piparav']/CR['niav']
        CR['tiperp1av'] = CR['piperp1av']/CR['niav']
        CR['tiperp2av'] = CR['piperp2av']/CR['niav']
        CR['teparav'] = CR['peparav']/CR['neav']
        CR['teperp1av'] = CR['peperp1av']/CR['neav']
        CR['teperp2av'] = CR['peperp2av']/CR['neav']

        print 'TOTAL TIME: %f'%(time.time() - t)

        CRL = {}
        CRU = {}
        ylen = len(CR['yy'])
        for k in CR:
            if k != 'yy' and k != 'xx':
                CRL[k] = np.squeeze(CR[k])[:ylen/2,:]
                CRU[k] = np.squeeze(CR[k])[ylen/2:,:]
            elif k == 'yy':
                CRL[k] = CR[k][:ylen/2]
                CRU[k] = CR[k][ylen/2:]
            else:
                CRL[k] = CR[k]
                CRU[k] = CR[k]
        
        print 'Saving lower data...'
        np.save(fname+'_lower',CRL)
        print 'Saving upper data...'
        np.save(fname+'_upper',CRU)

        if save_full:
            np.save(fname,CR)

    return CR


#======================================================

def roll_run(CR,sx=None):
    """ Roll every variable in a simulation (CR)
        in the x direction by length in indexspace (sx)
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
        if CR['yy'][0] < 1.0: 
            print 'Gonna roll RIGHT!!!!!!!!!!!!'
            sx = -1*np.size(CR['xx'])/4
        else: 
            print 'Gonna roll LEFT!!!!!!!!!!!!'
            sx = np.size(CR['xx'])/4
    
    for key in CR.keys():
# Old way not super smart
#        if key.rfind('av') == len(key)-2 and len(key) > 2:
#            print 'Rolling ',key
#            CR[key] = np.roll(CR[key],sx,axis=1)
# New way a little bit smarter
        if key in klst or key in kavlst:
            print 'Rolling ',key
            CR[key] = np.roll(CR[key],sx,axis=1)

#======================================================

def readsave(restore_fname):
    if restore_fname[restore_fname.rfind('.'):] == '.npy':
        cr = np.load(restore_fname).all()
        for v in ['ni', 'ne', 'niav', 'neav']:
            if v in cr:
                cr['de'+v] = cr[v]
        return cr
    else:
        return readsav(restore_fname)

#======================================================


def calc_pdf(ar,min=99999,max=99999,weight=100,inc=0,ax=0):
   if len(ar) == 0:
      print 'No array provided! Exiting!'
      return
   if min == 99999:
      min=ar.min()
   if max == 99999:
      max=ar.max()
   # If PDF of increment, then increment the array
   if inc > 0:
      ar = ar - np.roll(ar,inc,axis=ax)
   # Find the total length of data set
   arsize=reduce(operator.mul, np.shape(ar),1)
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

def rs3d(arr):
    """ Reshape an array as a 3D array (for Tulasi's stupid code)"""
    return arr.reshape(arr.shape + (1,))

def date_file_prefix():
    import datetime
    return datetime.date.today().strftime('%y.%m.%d')

#======================================================

def rotate_ten(CR,
               var='pi',
               av='av',
               overwrite=False,
               full_rotate=False):

    if var+'par'+av in CR and not overwrite:
        print 'Warning: %sparav was found in the' \
              'restored data: nothing will be rotated!!!!'
        pass

        
    elif full_rotate:
# e1 -> \hat{B} 
# e2 -> \hat{ExB}
# e3 -> \hat{Bx(ExB)}
        e1 = np.array([CR['bxav'],
                       CR['byav'],
                       CR['bzav']])

        e2 = np.cross(np.array([CR['exav'],
                                CR['eyav'],
                                CR['ezav']]), e1,axis=0)

        e1 = e1/np.sqrt(np.sum(e1**2,axis=0))
        e2 = e2/np.sqrt(np.sum(e2**2,axis=0))
        e3 = np.cross(e1,e2,axis=0)

        T = np.array([[CR[var+'xx'+av],CR[var+'xy'+av],CR[var+'xz'+av]],
                      [CR[var+'xy'+av],CR[var+'yy'+av],CR[var+'yz'+av]],
                      [CR[var+'xz'+av],CR[var+'yz'+av],CR[var+'zz'+av]]])

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
        CR[var+'par'+av]   = Tpar
        CR[var+'perp1'+av] = Tperp
        CR[var+'perp2'+av] = Tperp

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

        CR[var+'11'+av] = T11
        CR[var+'12'+av] = T12
        CR[var+'13'+av] = T13
        CR[var+'22'+av] = T22
        CR[var+'23'+av] = T23
        CR[var+'33'+av] = T23

        # Now Agyrotropy code!
        CR[var+'agy'+av] = 2.*np.sqrt((T12**2 + T13**2 + T23**2)) \
                                 /(T11 + T22 + T33) 
    else:
# This was the old way, and it was very simple

        bmag = np.sqrt( CR['bx'+av]**2+
                        CR['by'+av]**2+
                        CR['bz'+av]**2)
        bbx = CR['bx'+av]/bmag
        bby = CR['by'+av]/bmag
        bbz = CR['bz'+av]/bmag

        CR[var+'par'+av] = (bbx*(bbx*CR[var+'xx'+av] + 
                                 bby*CR[var+'xy'+av] + 
                                 bbz*CR[var+'xz'+av])+
                            bby*(bbx*CR[var+'xy'+av] + 
                                 bby*CR[var+'yy'+av] + 
                                 bbz*CR[var+'yz'+av])+
                            bbz*(bbx*CR[var+'xz'+av] + 
                                 bby*CR[var+'yz'+av] + 
                                 bbz*CR[var+'zz'+av]))

        CR[var+'perp1'+av] = (CR[var+'xx'+av] + 
                              CR[var+'yy'+av] +
                              CR[var+'zz'+av] - 
                              CR[var+'par'+av])/2.

        CR[var+'perp2'+av] = CR[var+'perp1'+av]

def date_file_prefix():
    import datetime
    return datetime.date.today().strftime('%Y.%m.%d.')
