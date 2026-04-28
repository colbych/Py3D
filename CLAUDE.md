# CLAUDE.md — Py3D Codebase Guide

## Project Overview

Py3D is a Python library for reading, analyzing, and visualizing particle-in-cell (PIC) plasma simulation data from the P3D simulation code. It is a research/academic codebase developed by Colby Haggerty, primarily targeting interactive use in IPython/Jupyter notebooks on HPC systems (e.g., NCAR Cheyenne/Yellowstone).

The codebase has been modernized (Python 3, pip-installable, pytest suite, ruff linting, non-interactive API) while keeping the research-library character intact. See `PLAN.md` for the active development roadmap.

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
├── tests/                   # pytest test suite
│   ├── conftest.py          # Shared fixtures
│   ├── test_methods.py      # Tests for py3d._methods
│   ├── test_sub.py          # Tests for py3d.sub utility functions
│   ├── test_vdist.py        # Tests for py3d.vdist.VDist
│   ├── test_movie.py        # Tests for py3d.movie.Movie (Phase 5)
│   └── test_dump.py         # Tests for py3d.dump.Dump and dumpID.DumpID (Phase 5)
│
├── pyproject.toml           # Package metadata, dependencies, tool config
├── README.txt               # Setup instructions
├── PLAN.md                  # Active development roadmap
└── CLAUDE.md                # This file
```

---

## Core Classes

### `py3d.Movie`
Loads time-series ("movie") simulation output.
- Supports two naming conventions: `'p3d'` (default) and `'tulasi'`
- Key method: `get_fields(t, *field_names)` — returns dict of numpy arrays at timestep `t`
- Auto-discovers parameter files via `load_param()`
- `interactive=True` (default) preserves original notebook behavior; pass `interactive=False` to raise instead of prompting on stdin (required for scripts/batch jobs)

### `py3d.Dump` (`dump.py`)
Low-level reader for P3D binary dump files.
- Parses custom binary format using Python's `struct` module
- Methods: `read_particles()`, `read_fields()`
- Handles both particle phase-space data and field data
- Accepts `interactive=False` for non-interactive use (same pattern as `Movie`)

### `py3d.DumpID` (`dumpID.py`)
Higher-level wrapper around `Dump` with spatial selection.
- Method `get_part_in_box(box)` — filters particles within a spatial box
- Supports coordinate rotation and transformation
- Accepts `interactive=False`, forwarded to the inner `Dump`

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
| `scipy` | Interpolation, IDL `.sav` file reading (`scipy.io.readsav`) |
| `matplotlib` | Plotting and visualization |
| `ctypes` | Interface to C extensions in `PartTrace/` |
| `mpi4py` | MPI parallelism (only in `DumpPartCompare/`) — optional `[hpc]` extra |

Dependencies are declared in `pyproject.toml`. Install with `pip install -e .`.

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
Interactive use (notebook/REPL — prompts for missing paths):
```python
import py3d
m = py3d.Movie('/path/to/simulation/')
d = m.get_fields(10, 'bx', 'by', 'bz', 'ex')
py3d.ims(d, 'bx')
```

Non-interactive use (scripts, batch jobs — raises instead of prompting):
```python
m = py3d.Movie(path='/path/to/sim/', num=0, param_file='p3d.in', interactive=False)
d = m.get_fields('bx by bz', time=10)
```

### Linting
```bash
ruff check py3d/ PartTrace/testparticle.py
```

Run this before committing. Configuration lives in `[tool.ruff]` in `pyproject.toml`. The rule set is `E` + `F` with a short ignore list for intentional style patterns (see the file for details).

### Running Tests
```bash
pip install -e ".[dev]"
pytest
```

Tests live in `tests/`. Pure functions and synthetic-data tests are fully covered. See **Deferred Testing Work** below for known gaps.

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

1. **Interactive prompts** *(mitigated in Phase 5)*: `Movie`, `Dump`, and `DumpID` still default to `interactive=True`, which calls `input()` when files are not found. Pass `interactive=False` to get a clean `FileNotFoundError`/`ValueError` instead. The `input()` fallback is preserved for backward compatibility in notebooks.
2. **C extensions**: `PartTrace` will silently fail or produce wrong results if `.so` files are not compiled for the current platform.
3. **Commented-out print statements**: Several files contain old Python 2 `print 'foo'` statements in comments — these are harmless but noisy. Do not mistake them for active code.

---

## Deferred Testing Work

Items still out of scope for the current test suite:

| Item | Reason deferred | Target phase |
|------|----------------|--------------|
| Plotting functions (`ims`, `PatPlotter.make_plots`, `VDistPlotter.plot2d`) | Produce matplotlib figures. Require visual regression tooling (e.g. `pytest-mpl`) to test meaningfully. | Post-Phase 5 |
| `PartTrace.TPRun` integration tests | Requires compiled C extensions (`functions.so`, `functions_rel.so`). Need a CI build step to compile them. | Phase 6 |
| `DumpPartCompare/` MPI tests | Requires an MPI runtime. Integration-test territory, not unit tests. | Out of scope |

Items completed in Phase 5:

| Item | Resolution |
|------|-----------|
| `Movie`, `Dump`, `DumpID` unit tests | Added `interactive=False` API; error-path and constructor tests in `tests/test_movie.py` and `tests/test_dump.py` |
| `VDist.spec1d` (full test) | Fixed `normed=True` → `density=True`; xfail marks removed; tests pass on NumPy 2.x |
