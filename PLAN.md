# Py3D — Development Plan

## What has been done

| Phase | Description | Status |
|-------|-------------|--------|
| 1 | Python 3 compatibility fixes | Done |
| 2 | pip-installable package (`pyproject.toml`) | Done |
| 3 | pytest test suite for pure functions | Done |
| 4 | ruff linting + bug fixes caught by static analysis | Done |

All work so far is on the open PR: https://github.com/colbych/Py3D/pull/15

---

## What is next: Phase 5 (split into steps)

The goal of Phase 5 is to make `Movie`, `Dump`, and `DumpID` testable and
usable non-interactively. Currently they call `input()` to ask for file paths
whenever a path is not found, which breaks use in scripts and batch jobs.

**Key constraint**: changes must be backward-compatible. Code that uses
`Movie()` today without arguments must continue to work.

### Step 5a — Fix `spec1d` (small, isolated, no risk)

`VDist.spec1d` calls `np.histogram2d(normed=True)`, which was removed in
NumPy 2.0. Change to `density=True`. Then remove the `xfail` marks from the
`test_spec1d_*` tests in `tests/test_vdist.py`.

File: `py3d/vdist.py` around line 202.

### Step 5b — Add optional path/num arguments to Movie, Dump, DumpID

Add keyword arguments so these classes can be instantiated without interactive
prompts, while keeping the `input()` fallback for existing interactive use.

Target API:
```python
# New: works non-interactively
m = Movie(path='/path/to/sim/', num=5, param_file='p3d.in')
d = Dump(path='/path/to/sim/', num=3)

# Old: still works, falls back to input() prompt as before
m = Movie()
```

Affected files and the specific `input()` call sites:
- `py3d/movie.py`: `_set_movie_path()` around line 324, `set_movie_num()` around line 346
- `py3d/dump.py`: `_set_dump_path()` around line 228, `set_dump_num()` around line 66
- `py3d/_methods.py`: `_get_param_file()` around line 47 (called by `load_param`)

The change pattern for each is the same: check if the argument was passed
before calling `input()`.

### Step 5c — Write tests for Movie, Dump, DumpID using the new API

Once Step 5b is done, write tests in `tests/test_movie.py` and
`tests/test_dump.py` using `tmp_path` fixtures and small synthetic binary
files (or minimal real param files).

The skipped stubs in `tests/test_methods.py` are the placeholder for this.

### Step 5d — Validate against real simulation data

**Before running Step 5b/5c code in production**, test the refactored classes
against actual P3D simulation data. The goal is to confirm that field arrays,
particle counts, and energy values from the new code path match the old
interactive path exactly.

Suggested validation script (run manually, not in pytest):
```python
import py3d

# Load with new explicit-path API
m = py3d.Movie(path='/path/to/sim/', num=0, param_file='p3d.in')
d = m.get_fields('bx by bz', time=10)

# Spot-check: print shapes, min/max of a few fields
for k, v in d.items():
    if hasattr(v, 'shape'):
        print(f"{k}: shape={v.shape}, min={v.min():.4f}, max={v.max():.4f}")
```

Run this before and after the refactor and compare outputs.

### Step 5e — Remove `input()` fallback (optional, later)

Once Step 5d passes, the `input()` fallback can be deprecated and eventually
removed. This is lower priority — the backward-compatible API from Step 5b is
good enough for most use cases.

---

## Deferred / out of scope

| Item | Reason | When to revisit |
|------|--------|----------------|
| Plotting tests (`ims`, `PatPlotter`, `VDistPlotter`) | Require visual regression tooling (`pytest-mpl`) | After Phase 5 |
| `PartTrace.TPRun` integration tests | Require compiled `.so` files and a CI build step | Phase 5 or later |
| `DumpPartCompare/` MPI tests | Require an MPI runtime; integration territory | Out of scope |
| Type hints | No annotations exist; adding them is a large effort for limited gain | Optional, low priority |
| f-string migration | Cosmetic; low value vs. effort on a working codebase | Optional |

---

## Notes for picking up this work

- Run `ruff check py3d/ PartTrace/testparticle.py` and `pytest` before and
  after any change to confirm nothing is broken.
- The binary file format for `Dump` is underdocumented. The struct layout is
  in `dump.py` around lines 479–492. Do not change parsing logic without
  validating against real data (Step 5d).
- `Movie` supports two naming conventions (`'p3d'` and `'tulasi'`). Any
  changes to path resolution need to be tested with both.
