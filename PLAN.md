# Py3D — Development Plan

## What has been done

| Phase | Description | Status |
|-------|-------------|--------|
| 1 | Python 3 compatibility fixes | Done |
| 2 | pip-installable package (`pyproject.toml`) | Done |
| 3 | pytest test suite for pure functions | Done |
| 4 | ruff linting + bug fixes caught by static analysis | Done |
| 5a | Fix `spec1d` NumPy 2.0 incompatibility (`normed` → `density`) | Done |
| 5b | Add `interactive=False` to `Movie`, `Dump`, `DumpID`, `load_param` | Done |
| 5c | Tests for `Movie`, `Dump`, `DumpID` using synthetic scaffolding | Done |
| 5d | Validate against real P3D simulation data (`py3d-validation/`) | Done |

Phases 1–4 merged via PR #15. Phase 5 work is on `master` (applied directly).

---

## What is next: Phase 6

Phase 5 is complete. Possible Phase 6 directions (priority TBD):

### Option A — `PartTrace.TPRun` CI build + tests
Add a `gcc` compile step for `PartTrace/functions.so` and `functions_rel.so`
so `TPRun` integration tests can run in CI. Currently skipped in the test
suite (`test_tprun_stub`).

### Option B — Plotting test coverage
Add `pytest-mpl` (or a simpler pixel-diff approach) to cover `ims`,
`PatPlotter.make_plots`, and `VDistPlotter.plot2d`. Low priority for a
research library but useful for catching regressions.

### Option C — Remove `input()` fallback (Step 5e, deferred)
Now that `interactive=False` is the recommended path, the `input()` loops
can be deprecated and eventually removed. Low urgency — the fallback is
harmless for interactive notebook use.

---

## Deferred / out of scope

| Item | Reason | When to revisit |
|------|--------|----------------|
| Plotting tests (`ims`, `PatPlotter`, `VDistPlotter`) | Require visual regression tooling (`pytest-mpl`) | Phase 6 Option B |
| `PartTrace.TPRun` integration tests | Require compiled `.so` files and a CI build step | Phase 6 Option A |
| `DumpPartCompare/` MPI tests | Require an MPI runtime; integration territory | Out of scope |
| Remove `input()` fallback | Low urgency; backward-compat is good enough | Phase 6 Option C |
| Type hints | No annotations exist; large effort for limited gain in a research lib | Optional |
| f-string migration | Cosmetic; low value vs. effort on a working codebase | Optional |

---

## Notes for picking up this work

- Run `ruff check py3d/ PartTrace/testparticle.py` and `pytest` before and
  after any change to confirm nothing is broken.
- Use `reconn-wave-power` conda env (NumPy 2.4, pytest 9) with
  `PYTHONPATH=/path/to/Py3D` for running tests without `pip install -e`.
- The binary file format for `Dump` is underdocumented. The struct layout is
  in `dump.py` around lines 479–492. Do not change parsing logic without
  validating against real data in `py3d-validation/`.
- `Movie` supports two naming conventions (`'p3d'` and `'tulasi'`). Any
  changes to path resolution need to be tested with both.
- Real-data validation scripts live in `~/Research/Programing/py3d-validation/`
  (outside the repo). Run `bash run_all.sh` there to validate against the
  example simulation at `P3D_example_simulation_data/staging/`.
