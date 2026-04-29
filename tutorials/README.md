# Py3D tutorials

Beginner-friendly examples for new users of Py3D, built around an
example asymmetric magnetic-reconnection simulation from the P3D code.

## Setup

See [SETUP.md](SETUP.md) for full instructions. Short version:

```bash
cd tutorials
conda env create -f environment.yml
conda activate py3d-tutorial
pip install -e ..
```

If the example dataset is not at the default path, set
`PY3D_TUTORIAL_DATA` to point at it.

## Tutorials

Each tutorial exists as both a runnable script and a hand-written
notebook. The two are structurally parallel — the notebook adds prose
and inline figures; the script is what the test suite runs.

| # | Topic | What you'll learn | Script | Notebook |
|---|-------|-------------------|--------|----------|
| 01 | Getting started | Load a movie, plot a field, overlay flux contours | [scripts/01_getting_started.py](scripts/01_getting_started.py) | [notebooks/01_getting_started.ipynb](notebooks/01_getting_started.ipynb) |
| 02 | Exploring fields | Load all 30 vars at one time, multi-panel plots | scripts/02_exploring_fields.py | notebooks/02_exploring_fields.ipynb |
| 03 | Time evolution | Loop over frames, build a reconnection-rate time series | scripts/03_time_evolution.py | notebooks/03_time_evolution.ipynb |
| 04 | Pressure tensor | Rotate into field-aligned frame, plot T∥ vs T⊥ | scripts/04_pressure_tensor.py | notebooks/04_pressure_tensor.ipynb |

Tutorials 05–06 (particle / VDist) are planned for a follow-up release.

## The example dataset

A 2D asymmetric reconnection run with these key parameters:

| Parameter | Value | Meaning |
|-----------|-------|---------|
| `lx`, `ly` | 51.2, 25.6 di | Domain size |
| `dt` | 0.01 Ωci⁻¹ | Timestep |
| `m_e` | 0.04 | mi/me = 25 |
| `b1`, `b2` | 1.0, 0.4 | Asymmetric upstream B |
| `n1`, `n2` | 1.0, 1.25 | Asymmetric upstream density |
| Frames | 16 | Movie snapshots available |

## Three ways to run a tutorial

Each script works equally well as a standalone program, an IPython
`%run` target, or a copy-paste source for a REPL.

**Standalone:**

```bash
python tutorials/scripts/01_getting_started.py
```

Runs top-to-bottom and saves a PNG next to the script.

**IPython interactive:**

```text
ipython --matplotlib
In [1]: %run tutorials/scripts/01_getting_started.py
In [2]: m            # Movie object — still in the namespace
In [3]: d.keys()     # the field dict from the script
In [4]: py3d.ims(d, 'by')   # poke at it some more
```

The scripts are flat (no `def main()` wrapper), so `%run` leaves every
top-level variable in your namespace for further exploration.

**Copy-paste / cell-by-cell:**

Each script is divided into sections by `# %%` markers. Paste blocks
into IPython one at a time, or use VS Code's "Run Cell" feature (the
Python extension recognises `# %%` automatically).

**Notebook:**

```bash
jupyter notebook tutorials/notebooks/01_getting_started.ipynb
```

## Conventions used in the tutorials

- **Non-interactive API everywhere.** All `Movie` and `Dump` instances
  are constructed with `interactive=False`, so missing files raise
  immediately instead of prompting on stdin.
- **Explicit paths.** `path` and `param_file` are always passed
  explicitly. `time=` is always passed to `get_fields`.
- **One data path source.** Scripts and notebooks both read `DATA_DIR`,
  `PARAM_FILE`, and `NAME_STYLE` from `data_path.py`.
