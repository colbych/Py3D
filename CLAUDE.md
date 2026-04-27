# CLAUDE.md — Py3D Codebase Guide

## Project Overview

Py3D is a Python library for reading, analyzing, and visualizing particle-in-cell (PIC) plasma simulation data from the P3D simulation code. It is a research/academic codebase developed by Colby Haggerty, primarily targeting interactive use in IPython/Jupyter notebooks on HPC systems (e.g., NCAR Cheyenne/Yellowstone).

There is no formal packaging, CI/CD, or test suite — this is intentional for a research library.

---

## Repository Structure

```
Py3D/
├── py3d/                    # Main source package
│   ├── __init__.py          # Exports top-level API
│   ├── dump.py              # Dump class — reads raw P3D dump files
│   ├── dumpID.py            # DumpID class — spatial filtering of dump data
│   ├── movie.py             # Movie class — time-series simulation data
│   ├── vdist.py             # VDist class — velocity distribution calculations
│   ├── vdist_plotter.py     # VDistPlotter class — velocity distribution plots
│   ├── patplots.py          # PatPlotter class — multi-page analysis plots
│   ├── sub.py               # Core utility functions (exported via *)
│   └── _methods.py          # Internal helpers (load_param, interp_field, etc.)
│
├── PartTrace/               # Test particle tracing module
│   ├── testparticle.py      # TPRun class — integrates particle trajectories
│   ├── functions.c          # C extension for fast field interpolation
│   └── functions_rel.c      # C extension for relativistic particle tracing
│
├── DumpPartCompare/         # MPI utility scripts
│   ├── part_dump_id.py      # Parallel particle-dump comparison (mpi4py)
│   ├── exe_partID.LSF       # LSF job submission script
│   └── exe_partID.cheyenne  # Cheyenne job submission script
│
├── README.txt               # Legacy setup instructions
└── CLAUDE.md                # This file
```

---

## Core Classes

### `py3d.Movie`
Loads time-series ("movie") simulation output.
- Supports two naming conventions: `'p3d'` (default) and `'tulasi'`
- Key method: `get_fields(t, *field_names)` — returns dict of numpy arrays at timestep `t`
- Auto-discovers parameter files via `load_param()`

### `py3d.Dump` (`dump.py`)
Low-level reader for P3D binary dump files.
- Parses custom binary format using Python's `struct` module
- Methods: `read_particles()`, `read_fields()`
- Handles both particle phase-space data and field data

### `py3d.DumpID` (`dumpID.py`)
Higher-level wrapper around `Dump` with spatial selection.
- Method `get_part_in_box(box)` — filters particles within a spatial box
- Supports coordinate rotation and transformation

### `py3d.VDist` / `py3d.VDistPlotter`
Velocity distribution analysis and plotting.
- `VDistPlotter` inherits from `VDist`
- Key methods: `vdist2d()`, `vdist2d_pitch()`, `spec1d()`, `plot2d()`

### `py3d.PatPlotter`
Multi-page figure generation for simulation analysis.
- Method `make_plots()` generates paginated analysis figures

### `PartTrace.TPRun`
Test particle trajectory integration.
- Uses ctypes to call compiled C extensions (`functions.c`, `functions_rel.c`) for performance
- Supports both relativistic and non-relativistic modes

---

## Key Utility Functions (`py3d.sub`)

Exported via `from py3d import *`:

| Function | Description |
|---|---|
| `ims(d, var)` | `imshow` wrapper for 2D simulation fields |
| `ims3D(d, var)` | 3D slice visualization |
| `load_movie(path)` | Convenience loader returning a `Movie` object |
| `load_parts(path)` | Load particle data |
| `calc_psi(d)` | Calculate magnetic flux function ψ |
| `find_xpt(d)` | Find X-points in reconnection simulations |
| `var_at(d, var, pos)` | Interpolate variable at a position |
| `multi_color(...)` | Multi-color line plotting by time index |
| `get_times(path)` | List available timesteps |
| `findval(arr, val)` | Index of 1D array element closest to `val` |
| `rotate_ten(...)` | Rotate a tensor |
| `show_energy(d)` | Display energy diagnostics |
| `calc_pdf(data)` | Calculate probability density function |

---

## Internal Helpers (`py3d._methods`)

| Function | Description |
|---|---|
| `load_param(path)` | Parse P3D `.in` parameter files into a dict |
| `vprint(msg, verbose)` | Verbose-conditional print |
| `interp_field(field, pos, d)` | Bilinear/trilinear field interpolation |
| `_num_to_ext(n)` | Convert integer to zero-padded 3-digit string |
| `_convert(val)` | String-to-int/float type coercion |
| `_get_param_file(path)` | Interactive discovery of parameter files |

---

## Dependencies

| Package | Purpose |
|---|---|
| `numpy` | Core numerical arrays |
| `scipy` | Interpolation, IDL `.sav` file reading (`scipy.io.idl`) |
| `matplotlib` | Plotting and visualization |
| `ctypes` | Interface to C extensions in `PartTrace/` |
| `mpi4py` | MPI parallelism (only in `DumpPartCompare/`) |

No `requirements.txt` exists. Install dependencies manually or via conda/pip.

---

## Code Conventions

### Style
- **Indentation**: 4 spaces
- **Class names**: PascalCase (`Movie`, `VDist`, `DumpID`)
- **Function/method names**: snake_case (`load_param`, `read_fields`)
- **Private/internal**: Leading underscore (`_methods.py`, `_read_fields_all`)
- **Short physics variables** are idiomatic: `r0`, `dx`, `v1`, `v2`, `v3`, `bx`, `by`, `bz`
- No type annotations — this is a pre-type-hints codebase; do not add them unnecessarily

### Docstrings
- Informal mix of reST and NumPy docstring styles
- Prefer documenting `Parameters` blocks, especially for public functions
- Example format used in the codebase:
  ```python
  """
  Short description.

  Parameters
  ==========
      d : dict
          Simulation data dictionary with xx, yy keys.
      var : str
          Field variable name to plot.
  """
  ```

### File Headers
All source files use this header format — preserve it when editing:
```python
#######################################################################
#                                                                     #
#                  Python Progs :  filename.py                        #
#                  Aruthor      :  Colby Haggerty                     #
#                  Date         :  YYYY.MM.DD                         #
#                                                                     #
#######################################################################
```

### Section Dividers
Use `#======...` style dividers between logical sections within files.

---

## Development Workflow

### Setup
Install as an editable package (recommended for research use — source edits take effect immediately):
```bash
pip install -e /path/to/Py3D
```

With optional MPI support for `DumpPartCompare/`:
```bash
pip install -e "/path/to/Py3D[hpc]"
```

If you cannot use pip, the legacy approach still works:
```bash
export PYTHONPATH=/path/to/Py3D:$PYTHONPATH
```

### Building C Extensions (PartTrace)
The `PartTrace` module requires compiled C extensions. Compile manually:
```bash
cd PartTrace
gcc -O3 -shared -fPIC -o functions.so functions.c
gcc -O3 -shared -fPIC -o functions_rel.so functions_rel.c
```

### Running Code
This library is designed for interactive use:
```python
import py3d
m = py3d.Movie('/path/to/simulation/')
d = m.get_fields(10, 'bx', 'by', 'bz', 'ex')
py3d.ims(d, 'bx')
```

### No Test Suite
There are no automated tests. When modifying code:
- Test interactively against real simulation data if available
- Check that `py3d/__init__.py` imports still work after changes
- Verify numpy array shapes and dtypes remain consistent

### Git Branching
- Main branch: `master`
- Feature/fix branches: use descriptive names, merge via pull request
- Recent contributor branches: `mshay/master`

---

## HPC / Supercomputer Context

- Job scripts in `DumpPartCompare/` target LSF (Yellowstone) and PBS (Cheyenne)
- `mpi4py` parallel code in `DumpPartCompare/part_dump_id.py` is intended for cluster runs
- Simulation data files can be very large (binary, multi-GB); memory-efficient access patterns are important

---

## Known Issues / Gotchas

1. **Python 2→3 migration**: The codebase was originally Python 2. Some legacy patterns may remain (e.g., `print` statement style comments, `input()` usage).
2. **No packaging**: Cannot be installed via `pip install .` — must use `PYTHONPATH`.
3. **Interactive prompts**: Some functions call `input()` to ask the user for file paths; avoid in non-interactive contexts.
4. **scipy.io.idl**: Used for reading IDL `.sav` files; may require older scipy versions depending on the simulation data format.
5. **C extensions**: `PartTrace` will silently fail or produce wrong results if `.so` files are not compiled for the current platform.
