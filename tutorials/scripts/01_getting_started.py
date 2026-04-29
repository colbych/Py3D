#######################################################################
#                                                                     #
#                  Python Progs :  01_getting_started.py              #
#                  Aruthor      :  Colby Haggerty                     #
#                  Date         :  2026.04.28                         #
#                                                                     #
#######################################################################
"""
Tutorial 01 — Getting started with Py3D.

The simplest possible Py3D session:
    1. Open a Movie pointing at the example simulation.
    2. Read one timestep of the magnetic field.
    3. Plot Bz with magnetic-flux contours overlaid.

Three ways to use this file:

    # Run end-to-end, saves a PNG next to the script:
    python tutorials/scripts/01_getting_started.py

    # Interactive: variables stay in your namespace afterwards:
    ipython --matplotlib
    In [1]: %run tutorials/scripts/01_getting_started.py
    In [2]: m            # the Movie object
    In [3]: d.keys()     # the field dict

    # Or paste each `# %%` block into IPython one at a time.
"""

# %% imports
import os
import sys

import matplotlib
import matplotlib.pyplot as plt

import py3d
from py3d.sub import calc_psi

# Add tutorials/ to sys.path so `from data_path import ...` works whether
# the script is run from the repo root or from tutorials/scripts/.
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from data_path import DATA_DIR, PARAM_FILE, NAME_STYLE, require_data_dir   # noqa: E402

require_data_dir()

# %% open the movie
# `interactive=False` makes Py3D raise on missing files instead of prompting
# on stdin. That's what you want for scripts. In a Jupyter session you can
# drop it; the default `interactive=True` will prompt for any path it
# cannot find.
m = py3d.Movie(
    num=0,                       # this run has only num=0
    path=DATA_DIR,
    param_file=PARAM_FILE,
    name_style=NAME_STYLE,
    interactive=False,
)
print(f'Loaded movie with {m.ntimes} frames')
print(f'Domain: lx = {m.param["lx"]}, ly = {m.param["ly"]}')

# %% read fields at one timestep
# `time` is a 0-based integer index into the available frames. Variable
# names go in as a single space-separated string. We grab bx, by, bz
# together so calc_psi has bx and by available later.
time_index = m.ntimes // 2     # midway through the run
d = m.get_fields('bx by bz', time=time_index)

print(f'Field shape: bz = {d["bz"].shape}')
print(f'             xx = {d["xx"].shape}, yy = {d["yy"].shape}')

# %% plot Bz with flux contours
# `py3d.ims` is the recommended way to display 2D fields. It transposes
# for you (the first axis of d['bz'] is x, not y), sets extents from
# xx/yy, and overlays psi contours when cont=True.
fig, ax = plt.subplots(figsize=(8, 4))
py3d.ims(d, 'bz', ax=ax, cbar=True, cont=True, cmap='RdBu_r')
ax.set_title(f'Bz at frame {time_index} (t = {d["tt"][-1]:.2f} Ωci⁻¹)')
fig.tight_layout()

out_png = os.path.join(HERE, '01_getting_started.png')
fig.savefig(out_png, dpi=120)
print(f'Saved {out_png}')

# %% reconnection diagnostic
# calc_psi returns the magnetic-flux function as a 2D array. The
# peak-to-peak excursion of psi is a quick proxy for how much flux has
# reconnected by this frame. Tutorial 03 plots this vs. time.
psi = calc_psi(d)
print(f'psi range at frame {time_index}: '
      f'[{psi.min():+.4f}, {psi.max():+.4f}], '
      f'span = {psi.max() - psi.min():.4f}')

# %% show the plot if running interactively
if matplotlib.get_backend().lower() not in ('agg', 'pdf', 'svg', 'ps'):
    plt.show()
