This anlysis package uses SciPy, NumPy and MatPlotLib.
If you are unfamilar with these, it might be good to go and read 
about them. Currently I am using this on yellowstone/cheyenne, hopefully 
nothing will change between super computers. So to use python on
yellowstone first you need to load the modules with;

module load python
module load all-python-libs

This will load all the nessesarty libraries
Then you can enter the interactive python front end
and run your scripts, with

ipython --pylab

You should also add the p3dthon direcotory to your path. I would just go ahead and add the following lines 
in your ipython config startup file:
in the file located at
    ~/.ipython/profile_default/startup/ipython_config.py 
Add the following code
#

import os
import sys

p3dthon_path = '/glade/u/home/colbyh/Py3D/'
sys.path.append(p3dthon_path)
from Py3D.sub import *
from Py3D.movie import Movie
from Py3D.vdist_plotter import VDistPlotter
from PartTrace.testparticle import TPRun

#
But replase the path with your own path to Py3D

this will let you use Py3D right off the bat.
Good luck!

