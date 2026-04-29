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
| 6a | Tutorials — env setup, scaffolding, 4 field/movie tutorials, smoke test | Done |

Phases 1–4 merged via PR #15. Phase 5 work is on `master` (applied directly).

---

## What is next: Phase 6 — Tutorials

Add a `tutorials/` directory with beginner-friendly examples for new users.
Driven by the example dataset at
`/Users/colby/Research/Programing/P3D_example_simulation_data/staging/`
(param file `param_asym00p`, `name_style='p3d'`, asymmetric reconnection,
16 movie frames, 16 dump snapshots).

Reference: `~/Research/Programing/py3d-validation/TUTORIAL_NOTES.md` —
captures the API patterns and gotchas confirmed against real data.

### Layout

```
tutorials/
├── README.md                # tutorial index + 4-line pointer to SETUP.md
├── SETUP.md                 # conda env walkthrough, py3d editable install, data path
├── environment.yml          # one-shot `conda env create -f` spec (loose pins)
├── data_path.py             # central DATA_DIR + PARAM_FILE; overridable via PY3D_TUTORIAL_DATA
├── scripts/                 # plain .py — runnable end-to-end, covered by smoke test
└── notebooks/               # hand-written .ipynb, structurally parallel to scripts
```

### Step 0 — Particle smoke test (prerequisite for 6b)  ✅ DONE

Throwaway script `~/Research/Programing/py3d-validation/validate_particles_staging.py`
confirms `Dump(...).read_particles(index=1)` works against the staging
data. Configuration: `pex*pey*pez = 256`, `nchannels = 16`, divisibility OK.
Quasi-neutrality holds (Ne/Ni = 0.9998), positions inside domain, vth_e/vth_i
ratio matches mi/me=25 expectation.

**Observation for Phase 6b**: each single-`num` dump file holds only a y-slab
of the domain (16-proc band of the 16×16 proc grid), not the full grid.
Tutorial 05 should use `DumpID.get_part_in_box` to aggregate across `num`
values — that's `DumpID`'s purpose anyway and a real reason for new users
to learn it.

### Phase 6a — Field/movie tutorials (one PR)  ✅ DONE

Four tutorials covering everything reachable through `Movie`:

| # | Title | API surface |
|---|-------|-------------|
| 01 | Getting started | `Movie`, `get_fields`, `ims`, `calc_psi` |
| 02 | Exploring fields | `get_fields('all', time=t)`, `movie_vars`, multi-axes `ims` |
| 03 | Time evolution | loop over `ntimes`, `find_xpt`, reconnection rate |
| 04 | Pressure tensor | `rotate_ten`, T∥/T⊥ maps |

Plus:
- `tutorials/README.md` (index + setup pointer) and `tutorials/data_path.py`
- `tutorials/SETUP.md` and `tutorials/environment.yml` — conda env walkthrough
- `tests/test_tutorials.py` — runs each script if the data dir exists,
  skips otherwise; marked `slow` so default `pytest` stays fast

### Phase 6b — Particle/VDist tutorials (separate PR)

Gated on the Step 0 smoke test passing.

| # | Title | API surface |
|---|-------|-------------|
| 05 | Velocity distributions | `Dump.read_particles`, `DumpID.get_part_in_box`, `VDist.vdist2d` |
| 06 | Energy spectrum | `VDist.spec1d` (NumPy-2.x compatible since Phase 5a) |

Extend `tests/test_tutorials.py` to cover 05–06.

### Recorded design decisions

- **Notebook format**: hand-written `.ipynb` (not jupytext-generated).
  Six small notebooks aren't worth a build step, and what's checked in is
  exactly what a new user opens in Jupyter. Scripts and notebooks are
  structurally parallel but not byte-synced.
- **Non-interactive API everywhere**: tutorials always use `interactive=False`
  and pass `path`, `param_file`, `time` explicitly. Failures are loud.
- **No new runtime deps**. Notebooks need `jupyter` to run, but that's the
  user's responsibility — noted in the README.
- **Scripts are flat**, not wrapped in `def main()` / `if __name__ == '__main__'`.
  This makes `%run scripts/0X_*.py` in IPython leave `m`, `d`, etc. in the
  user's namespace, and lets users copy-paste blocks directly into a REPL.
  `# %%` markers separate logical sections (also recognised by VS Code's
  Python interactive mode). Tutorials are not meant to be imported as
  modules, so the loss of the `main()` convention is a non-issue.

---

## Phase 7+ candidates (deferred, no commitment)

### `PartTrace.TPRun` CI build + tests
Add a `gcc` compile step for `PartTrace/functions.so` and `functions_rel.so`
so `TPRun` integration tests can run in CI. Currently skipped in the test
suite (`test_tprun_stub`).

### Plotting test coverage
Add `pytest-mpl` (or a simpler pixel-diff approach) to cover `ims`,
`PatPlotter.make_plots`, and `VDistPlotter.plot2d`. Low priority for a
research library but useful for catching regressions.

### Remove `input()` fallback
Now that `interactive=False` is the recommended path, the `input()` loops
can be deprecated and eventually removed. Low urgency — the fallback is
harmless for interactive notebook use.

---

## Deferred / out of scope

| Item | Reason | When to revisit |
|------|--------|----------------|
| Plotting tests (`ims`, `PatPlotter`, `VDistPlotter`) | Require visual regression tooling (`pytest-mpl`) | Phase 7+ |
| `PartTrace.TPRun` integration tests | Require compiled `.so` files and a CI build step | Phase 7+ |
| `DumpPartCompare/` MPI tests | Require an MPI runtime; integration territory | Out of scope |
| Remove `input()` fallback | Low urgency; backward-compat is good enough | Phase 7+ |
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
