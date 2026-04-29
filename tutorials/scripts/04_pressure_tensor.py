#######################################################################
#                                                                     #
#                  Python Progs :  04_pressure_tensor.py              #
#                  Aruthor      :  Colby Haggerty                     #
#                  Date         :  2026.04.28                         #
#                                                                     #
#######################################################################
"""
Tutorial 04 — Pressure tensor and temperature anisotropy.

Loads the full pressure tensor (six components per species) and uses
py3d.sub.rotate_ten to project it into the field-aligned frame. Then
computes the parallel and perpendicular temperatures and plots the
anisotropy T_par / T_perp.

Parallel/perpendicular here is "with respect to B", not the simulation
axes — that's what `rotate_ten` does for you.

See the header of 01_getting_started.py for the three ways to run this.
"""

# %% imports
import os
import sys

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import py3d
from py3d.sub import rotate_ten

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from data_path import DATA_DIR, PARAM_FILE, NAME_STYLE, require_data_dir   # noqa: E402

require_data_dir()

# %% open the movie
m = py3d.Movie(
    num=0,
    path=DATA_DIR,
    param_file=PARAM_FILE,
    name_style=NAME_STYLE,
    interactive=False,
)

# %% load fields needed for both species
# Each call to rotate_ten needs the magnetic field plus all six
# components of one pressure tensor. We grab densities too so we can
# divide pressure by density to get temperature.
time_index = m.ntimes // 2
ion_keys = 'pixx piyy pizz pixy piyz pixz ni'
ele_keys = 'pexx peyy pezz pexy peyz pexz ne'
b_keys = 'bx by bz'

d = m.get_fields(f'{b_keys} {ion_keys} {ele_keys}', time=time_index)
print(f'Loaded fields at frame {time_index}, shape = {d["bx"].shape}')

# %% rotate both pressure tensors into the field-aligned frame
# rotate_ten adds {var}par, {var}perp1, {var}perp2 to d in place.
# av='' tells it to use 'bx', 'by', 'bz' (no 'av' suffix). var='pi' or
# 'pe' selects the species.
rotate_ten(d, var='pi', av='')
rotate_ten(d, var='pe', av='')

print('After rotate_ten, new keys:')
for key in ('pipar', 'piperp1', 'pepar', 'peperp1'):
    arr = d[key]
    print(f'  {key:>8s}  range = [{arr.min():.3f}, {arr.max():.3f}]')

# %% compute parallel and perpendicular temperatures
# Temperature is pressure / density. piperp1 == piperp2 in the simple
# rotation (the axisymmetry assumption), so just use perp1.
Ti_par = d['pipar'] / d['ni']
Ti_perp = d['piperp1'] / d['ni']
Te_par = d['pepar'] / d['ne']
Te_perp = d['peperp1'] / d['ne']

# Anisotropy: T_par / T_perp. == 1 is isotropic;
# > 1 is firehose-prone, < 1 is mirror-prone.
ion_aniso = Ti_par / Ti_perp
ele_aniso = Te_par / Te_perp

print('\nAnisotropy (T_par / T_perp) ranges:')
print(f'  ions:      [{ion_aniso.min():.3f}, {ion_aniso.max():.3f}], '
      f'mean = {ion_aniso.mean():.3f}')
print(f'  electrons: [{ele_aniso.min():.3f}, {ele_aniso.max():.3f}], '
      f'mean = {ele_aniso.mean():.3f}')

# %% pack T_par, T_perp, anisotropy into the dict so ims can plot them
# `py3d.ims` reads the variable straight out of `d` by name.
d['Ti_par'] = Ti_par
d['Ti_perp'] = Ti_perp
d['ion_aniso'] = ion_aniso
d['Te_par'] = Te_par
d['Te_perp'] = Te_perp
d['ele_aniso'] = ele_aniso

# %% plot 2x3: T_par, T_perp, anisotropy for ions and electrons
fig, axes = plt.subplots(
    nrows=2, ncols=3, figsize=(15, 7), sharex=True, sharey=True
)

panels = [
    [('Ti_par', 'Ti par (ions)', None),
     ('Ti_perp', 'Ti perp (ions)', None),
     ('ion_aniso', 'Ti par / Ti perp', LogNorm(vmin=0.5, vmax=2.0))],
    [('Te_par', 'Te par (electrons)', None),
     ('Te_perp', 'Te perp (electrons)', None),
     ('ele_aniso', 'Te par / Te perp', LogNorm(vmin=0.5, vmax=2.0))],
]

for i, row in enumerate(panels):
    for j, (var, title, norm) in enumerate(row):
        ax = axes[i, j]
        if norm is not None:
            py3d.ims(d, var, ax=ax, cbar=True, cmap='RdBu_r',
                     norm=norm)
        else:
            py3d.ims(d, var, ax=ax, cbar=True, cmap='inferno')
        ax.set_title(title)
        if i == 0:
            ax.set_xlabel('')
        if j > 0:
            ax.set_ylabel('')

fig.suptitle(
    f'Field-aligned temperatures at frame {time_index} '
    f'(t = {time_index * m.param["n_movieout"] * m.param["dt"]:.2f} Ωci⁻¹)',
    fontsize=14,
)
fig.tight_layout()

out_png = os.path.join(HERE, '04_pressure_tensor.png')
fig.savefig(out_png, dpi=120)
print(f'Saved {out_png}')

# %% summary numbers
ion_iso_frac = float(np.mean(np.abs(np.log(ion_aniso)) < 0.1))
ele_iso_frac = float(np.mean(np.abs(np.log(ele_aniso)) < 0.1))
print('\nFraction of grid within 10% of isotropic '
      '(|log(T_par/T_perp)| < 0.1):')
print(f'  ions:      {ion_iso_frac:6.1%}')
print(f'  electrons: {ele_iso_frac:6.1%}')

# %% show the figure if running interactively
if matplotlib.get_backend().lower() not in ('agg', 'pdf', 'svg', 'ps'):
    plt.show()
