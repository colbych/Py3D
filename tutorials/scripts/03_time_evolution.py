#######################################################################
#                                                                     #
#                  Python Progs :  03_time_evolution.py               #
#                  Aruthor      :  Colby Haggerty                     #
#                  Date         :  2026.04.28                         #
#                                                                     #
#######################################################################
"""
Tutorial 03 — Time evolution and reconnection rate.

Loops over every available frame, computes two simple reconnection
diagnostics per frame, and plots them as a time series:

    - psi excursion  (max - min of the magnetic flux function)
    - peak |Ez|      (a proxy for the reconnection electric field)

Reading 51 frames takes ~1 minute on the example dataset. If you want
something faster, change FRAME_STRIDE.

See the header of 01_getting_started.py for the three ways to run this.
"""

# %% imports
import os
import sys

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import py3d
from py3d.sub import calc_psi

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

# Stride lets you skip frames if you don't want to wait for all of them.
# Set to 1 for the full set, 2 for every other frame, etc.
FRAME_STRIDE = 1
frames = list(range(0, m.ntimes, FRAME_STRIDE))
print(f'Will read {len(frames)} of {m.ntimes} frames')

# %% time axis from sim parameters
# Movie output cadence: every n_movieout simulation steps of length dt.
dt_movie = m.param['n_movieout'] * m.param['dt']
times = np.array(frames, dtype=float) * dt_movie

# %% loop over frames, compute diagnostics
# Reading several variables at once is faster than reading them one by
# one — the file is opened once per get_fields call.
psi_span = np.empty(len(frames))
ez_peak = np.empty(len(frames))

for i, t in enumerate(frames):
    df = m.get_fields('bx by ez', time=t)
    psi = calc_psi(df)
    psi_span[i] = psi.max() - psi.min()
    ez_peak[i] = np.abs(df['ez']).max()
    if (i + 1) % 10 == 0 or i == len(frames) - 1:
        print(f'  frame {t:>3d}  t = {times[i]:6.2f}   '
              f'psi span = {psi_span[i]:6.3f}   '
              f'peak |Ez| = {ez_peak[i]:6.3f}')

# %% reconnected-flux growth rate
# Reconnection rate is conventionally the time derivative of the flux
# transferred across the X-line. We approximate that here by the time
# derivative of the psi excursion, smoothed over neighbouring frames
# via numpy.gradient.
rec_rate = np.gradient(psi_span, times)

# %% two-panel time series
fig, (ax1, ax2) = plt.subplots(
    nrows=2, ncols=1, figsize=(8, 6), sharex=True
)

ax1.plot(times, psi_span, 'k-', lw=1.5, label='psi span')
ax1b = ax1.twinx()
ax1b.plot(times, rec_rate, 'C3--', lw=1.0, label='d(psi span)/dt')
ax1.set_ylabel('psi excursion (max − min)', color='k')
ax1b.set_ylabel('rate of change', color='C3')
ax1b.tick_params(axis='y', colors='C3')
ax1.set_title('Reconnected flux and its growth rate')

ax2.plot(times, ez_peak, 'C0-', lw=1.5)
ax2.set_xlabel('t  (Ωci⁻¹)')
ax2.set_ylabel('peak |Ez|')
ax2.set_title('Peak reconnection electric field')

for ax in (ax1, ax2):
    ax.grid(True, alpha=0.3)

fig.tight_layout()
out_png = os.path.join(HERE, '03_time_evolution.png')
fig.savefig(out_png, dpi=120)
print(f'Saved {out_png}')

# %% summary
peak_idx = int(np.argmax(rec_rate))
print(f'\nPeak reconnection rate {rec_rate[peak_idx]:.4f} '
      f'at t = {times[peak_idx]:.2f} (frame {frames[peak_idx]})')
print(f'Peak |Ez| over the run = {ez_peak.max():.4f} '
      f'at t = {times[int(np.argmax(ez_peak))]:.2f}')

# %% show the figure if running interactively
if matplotlib.get_backend().lower() not in ('agg', 'pdf', 'svg', 'ps'):
    plt.show()
