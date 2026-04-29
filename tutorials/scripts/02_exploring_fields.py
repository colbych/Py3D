#######################################################################
#                                                                     #
#                  Python Progs :  02_exploring_fields.py             #
#                  Aruthor      :  Colby Haggerty                     #
#                  Date         :  2026.04.28                         #
#                                                                     #
#######################################################################
"""
Tutorial 02 — Exploring all the fields.

Loads every variable that the simulation wrote out at a single timestep
and arranges nine of them (B, E, J vector components) in a 3x3 panel
plot. Along the way:

    - List m.movie_vars to see what's available.
    - Use get_fields('all', time=t) to read everything in one call.
    - Reuse the same axis layout across many ims calls.

See the header of 01_getting_started.py for the three ways to run this.
"""

# %% imports
import os
import sys

import matplotlib
import matplotlib.pyplot as plt

import py3d

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from data_path import DATA_DIR, PARAM_FILE, NAME_STYLE, require_data_dir   # noqa: E402

require_data_dir()

# %% open the movie and list available variables
m = py3d.Movie(
    num=0,
    path=DATA_DIR,
    param_file=PARAM_FILE,
    name_style=NAME_STYLE,
    interactive=False,
)
print(f'{len(m.movie_vars)} movie variables:')
print('  ' + ', '.join(m.movie_vars))

# %% read everything at one timestep
# Passing 'all' is shorthand for the full list. The return is a dict
# keyed by variable name plus the usual xx, yy, tt, time, param.
time_index = m.ntimes // 2
d = m.get_fields('all', time=time_index)

# Sanity check: every named variable made it in
missing = [v for v in m.movie_vars if v not in d]
assert not missing, f'missing from d: {missing}'
print(f'Read {len(m.movie_vars)} fields at frame {time_index} '
      f'(t = {d["tt"][-1]:.2f} Ωci⁻¹)')

# %% pick a 3x3 set of vector components
# Three rows: B, E, J — the dominant fields in MHD-scale reconnection.
# Three columns: x, y, z components.
panels = [
    ['bx', 'by', 'bz'],
    ['ex', 'ey', 'ez'],
    ['jx', 'jy', 'jz'],
]

# %% build the multi-panel figure
# The same `ims` call we used in tutorial 01, just looped. We pass cont=True
# only on the top-left panel so the flux contours appear once for context
# without cluttering every subplot.
fig, axes = plt.subplots(
    nrows=3, ncols=3, figsize=(14, 8), sharex=True, sharey=True
)

for i, row in enumerate(panels):
    for j, var in enumerate(row):
        ax = axes[i, j]
        py3d.ims(d, var, ax=ax, cbar=True,
                 cont=(i == 0 and j == 0),
                 cmap='RdBu_r')
        ax.set_title(var)
        # Only label outer axes to keep the panel clean.
        if i < 2:
            ax.set_xlabel('')
        if j > 0:
            ax.set_ylabel('')

fig.suptitle(f'B, E, J at frame {time_index} (t = {d["tt"][-1]:.2f} Ωci⁻¹)',
             fontsize=14)
fig.tight_layout()

out_png = os.path.join(HERE, '02_exploring_fields.png')
fig.savefig(out_png, dpi=120)
print(f'Saved {out_png}')

# %% quick comparison: peak amplitude per variable
# A useful one-liner once everything is loaded — sometimes the easiest
# way to know what's "going on" in a simulation is to see which fields
# have the most signal.
print('\nPeak |value| by variable:')
for var in m.movie_vars:
    arr = d[var]
    print(f'  {var:>5s}  max |.| = {abs(arr).max():.4f}   '
          f'std = {arr.std():.4f}')

# %% show the figure if running interactively
if matplotlib.get_backend().lower() not in ('agg', 'pdf', 'svg', 'ps'):
    plt.show()
