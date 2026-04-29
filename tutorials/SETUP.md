# Tutorial setup

This guide gets you a working environment for the Py3D tutorials. Plan
on five minutes the first time, plus however long conda takes to solve.

## 1. Prerequisites

You need a conda distribution. Any of these works:
[Miniconda](https://docs.conda.io/en/latest/miniconda.html),
[Mambaforge](https://github.com/conda-forge/miniforge),
or full Anaconda. The instructions below use `conda`; substitute `mamba`
if you have it (it is faster).

## 2. Create the environment

From the `tutorials/` directory:

```bash
cd /path/to/Py3D/tutorials
conda env create -f environment.yml
conda activate py3d-tutorial
```

This installs Python, NumPy, SciPy, matplotlib, and Jupyter. Pins are
loose — current versions are expected to work.

## 3. Install Py3D in editable mode

From the same shell (in the root of the py3d directory), 
with the env active:

```bash
pip install -e ..
```

Editable (`-e`) means edits to the Py3D source under `../py3d/` take
effect immediately the next time you import it. This is the recommended
mode for tutorial use, since you will probably want to read or tweak the
library code as you learn.

## 4. Verify the install

```bash
python -c "import py3d; print('py3d OK from', py3d.__file__)"
```

You should see a path ending in `Py3D/py3d/__init__.py`. Anything else
means the import picked up a different copy — check that the env is
active.

## 5. Point the tutorials at your data

The tutorials read from a fixed default path:

```
/Users/colby/Research/Programing/P3D_example_simulation_data/staging
```

If your example dataset lives somewhere else, set:

```bash
export PY3D_TUTORIAL_DATA=/path/to/your/staging
```

(Add this to your shell profile if you want it permanent.) The
`tutorials/data_path.py` module reads this env var and falls back to the
default path otherwise.

## 6. Run a tutorial

Scripts:

```bash
python scripts/01_getting_started.py
```

Notebooks:

```bash
jupyter notebook notebooks/01_getting_started.ipynb
```

## Troubleshooting

**`ModuleNotFoundError: No module named 'py3d'`**
The env is not active, or `pip install -e ..` was run from the wrong
directory. Re-run step 3 from inside `tutorials/`.

**`FileNotFoundError: Tutorial data directory not found`**
The default path does not exist on your machine. Set
`PY3D_TUTORIAL_DATA` (step 5).

**Plots do not appear**
For scripts: matplotlib falls back to the default backend, which on
some headless systems is non-interactive. Try `MPLBACKEND=TkAgg` or run
the notebook version instead. The scripts also write a PNG copy of each
figure next to the script, so you can view that even without a display.

**`conda env create` is slow**
Use `mamba env create -f environment.yml` instead of `conda env create`.
Mamba ships with Mambaforge or can be installed into any conda env via
`conda install -n base -c conda-forge mamba`.
