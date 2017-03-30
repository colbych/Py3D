import os
import numpy as np

def load_param(param_file=None):
    """ Method to load in the param file for a given run
        It will try and then ask for where the file is. if it doent know
    """
# Add a try catch statment incase you cant find the file

    param = {'file': param_file}

    if param['file'] is None:
        param['file'] = _get_param_file()

    fname = os.path.abspath(os.path.expandvars(param['file']))

    with open(fname) as f:
        content = f.readlines()

    for item in content:
        if '#define' in item and item[0] != '!':
            if len(item.split()) > 2:
                key = item.split()[1]
                val = item.split()[2]
                val = _convert(item.split()[2])
            else:
                key = item.split()[1]
                val = None

            param[key] = val

# An issue: In the param it is commen to say nchannels as pex. This
#           presents a problem in how we read the param file. So we
#           can run through the param dictionary and replace any value
#           with the coresponding key
    
    for key,val in param.iteritems():
        if param.has_key(val):
            param[key] = param[val]

    return param


def _get_param_file():
    fname = raw_input('Please Param File: ')
    fname = os.path.abspath(os.path.expandvars(fname.strip()))

    while not os.path.isfile(fname):
        error_text = '\nFile %s not found!\n' \
                     'Please Enter Param file: ' % fname

        fname = raw_input(error_text)
        fname = os.path.abspath(os.path.expandvars(fname))

    return fname


def _convert(val):
    constructors = [int, float, str]
    for c in constructors:
        try:
            return c(val)
        except ValueError:
            pass


def _num_to_ext(num):
    if num is not None:
        return '{0:03d}'.format(int(num))
    else:
        return None

###### This should be in a differnt class

def interp_field(fld, r0, sim_lens):
    r0 = np.array(r0)
    sim_lens = np.array(sim_lens)
    nl = np.array(np.shape(fld))
    
    dl = 1.0*sim_lens/nl

    lp = np.floor((r0 - dl/2.)/sim_lens*nl).astype(int)
    
    lp1 = (lp + 1)%nl

    wl = (r0 - (lp + .5)*dl)/dl

    for w in wl:
        if w < 0.:
            print 'C'*80
            print 'Negetive interp weight'
            print 'C'*80

    if len(np.shape(fld)) == 2:
        return (1.-wl[0])*(1.-wl[1])*fld[lp[0],  lp[1] ]+\
               (wl[0])   *(1.-wl[1])*fld[lp[0],  lp1[1]]+\
               (1.-wl[0])*   (wl[1])*fld[lp1[0], lp[1] ]+\
               (wl[0])   *   (wl[1])*fld[lp1[0], lp1[1]]

    elif len(np.shape(fld)) == 3:
        return (1.-wl[0])*(1.-wl[1])*(1.-wl[2])*fld[lp[0],  lp[1],  lp[2] ]+\
               (wl[0])   *(1.-wl[1])*(1.-wl[2])*fld[lp[0],  lp1[1], lp[2] ]+\
               (1.-wl[0])*   (wl[1])*(1.-wl[2])*fld[lp1[0], lp[1],  lp[2] ]+\
               (wl[0])   *   (wl[1])*(1.-wl[2])*fld[lp1[0], lp1[1], lp[2] ]+\
               (1.-wl[0])*(1.-wl[1])*   (wl[2])*fld[lp[0],  lp[1],  lp1[2]]+\
               (wl[0])   *(1.-wl[1])*   (wl[2])*fld[lp[0],  lp1[1], lp1[2]]+\
               (1.-wl[0])*   (wl[1])*   (wl[2])*fld[lp1[0], lp[1],  lp1[2]]+\
               (wl[0])   *   (wl[1])*   (wl[2])*fld[lp1[0], lp1[1], lp1[2]]
    else:
        raise Exception("Field shape not understood")
