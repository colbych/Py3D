Py3D — Python tools for P3D particle-in-cell simulation analysis
=================================================================

Dependencies: NumPy, SciPy, Matplotlib

Installation
------------
Install as an editable package (recommended — source edits take effect immediately):

    pip install -e /path/to/Py3D

With optional MPI support for DumpPartCompare/:

    pip install -e "/path/to/Py3D[hpc]"

On HPC systems without pip access, add the repo to your Python path instead:

    export PYTHONPATH=/path/to/Py3D:$PYTHONPATH

Basic usage
-----------
    import py3d
    m = py3d.Movie('/path/to/simulation/')
    d = m.get_fields(10, 'bx', 'by', 'bz', 'ex')
    py3d.ims(d, 'bx')

PartTrace (test particle tracing)
----------------------------------
PartTrace requires compiled C extensions. Build them manually before use:

    cd PartTrace
    gcc -O3 -shared -fPIC -o functions.so functions.c
    gcc -O3 -shared -fPIC -o functions_rel.so functions_rel.c

Then:

    from PartTrace.testparticle import TPRun

HPC modules (Cheyenne/Yellowstone)
------------------------------------
    module load python
    module load all-python-libs
