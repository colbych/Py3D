#!/usr/bin/env python
"""
Particle Dump File Compare

This code is run from the script exe_partID.LSF
In that script you need to pass the arguments as a colon ':'
sperated string where the even args are the var name and the odd args
will be the values.
   r0: (list) center of box (x,y,z) of particles we are collecting 
   dx: (list) width, height and depth of the box of particles we are 
       collecting 
   sp: ('i' or 'e') species of particles to collect and compare
   param_file: (str) name of the param file for the run
   init_path: (str) path of the dumpfiles you want to search through
   init_dump_num: (int) numeric extension for the dump files
   end_path: (str) path of the dumpfiles you want to collect in the box 
   end_dump_num: (int) numeric extension for the end dump files
"""
from mpi4py import MPI
import sys
import numpy as np

p3dthon_path = '/glade/u/home/colbyh/Py3D/'
sys.path.append(p3dthon_path)

from Py3D.dump import Dump
from Py3D.dumpID import DumpID

############ VALUES ARE SET IN THE EXE FILE! ###############
############  please don't edit this unless  ###############
############  you know what you are doing!   ###############

def cast_arg(v):
    """ A little method to convert the arguments from the exe file
    """
    if v[0] == '[' and v[-1] == ']':
        vals = v[1:-1].split(',')
        return [float(k) for k in vals]

    for tp in [int, str]:
        try:
            return tp(v)
        except ValueError:
            continue

arg = ''
for a in sys.argv[1:]:
    arg+=a

arg = arg.split(':')
for k,v in zip(arg[::2], arg[1::2]):
    locals()[k] = cast_arg(v)


save_name = 'init_partID_{}_r0=['+ (len(r0)*'{:.2f}, ')[:-2] + \
            ']_dx=['+(len(r0)*'{:.2f}, ')[:-2] + ']_{}'

save_name = save_name.format(*([param_file]+r0+dx+[sp]))

size = MPI.COMM_WORLD.Get_size()
rank = MPI.COMM_WORLD.Get_rank()

# Read the particles in the Box we care about
if rank == 0:
    did = DumpID(num=end_dump_num, param_file=param_file, path=end_path)
    parts = did.get_part_in_box(r0,dx,tags=True)
    data = parts[sp]
    param = did.param
    final_parts = None

else:
    data = None
    param = None

# broadcast that info to all of the processors 
param = MPI.COMM_WORLD.bcast(param,root=0)
data = MPI.COMM_WORLD.bcast(data,root=0)

if param['nchannels'] != size:
    print 'Compiled for wrong number of processors!!! Exiting!!!'
#    sys.exit()

D = Dump(num=init_dump_num, param_file=param_file, path=init_path)

sub_array_on_dump = 1
for k in 'xyz':
    sub_array_on_dump*= D.param['pe'+k]

sub_array_on_dump/=D.param['nchannels']
if rank == 0:
    print 'sub_array_on_dump ', sub_array_on_dump

MPI.COMM_WORLD.Barrier()

found_parts = []
for c in range(sub_array_on_dump):
    temp_parts = D.read_particles(rank+1,[c],tags=True)
    tps = temp_parts[sp][c]
    matches = np.in1d(tps['tag'], data['tag'], assume_unique=True)
    if True in matches:
        found_parts.append(tps[matches])

    if rank == 0:
        print 'PE-{}: Done {}/{}'.format(rank,c,sub_array_on_dump)

MPI.COMM_WORLD.Barrier()
if rank == 0:
    print 'Trying to gather the found particles'
final_parts = MPI.COMM_WORLD.gather(found_parts, root=0)
if rank == 0:
    fps = []
    for k in final_parts:
        fps += k
    print 'Saving found parts',rank+1 
    np.save(save_name, dict(found_parts=np.concatenate(fps),
                            parts=parts,r0=r0,dx=dx))
