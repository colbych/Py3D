"""
test_tutorials.py — smoke tests for the tutorial scripts.

Runs each tutorial script as a subprocess against the staging dataset
and confirms it exits cleanly. Skipped when the data directory is not
present, so this is a no-op on machines without the example data.

Marked `slow` because the longest scripts read all 51 movie frames
(a few seconds, but enough that we don't want it on the default
`pytest` run). Enable with:

    pytest -m slow

Or run just one script:

    pytest -m slow tests/test_tutorials.py::test_tutorial_script -k 03
"""

import os
import subprocess
import sys

import pytest


REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TUTORIAL_SCRIPTS = os.path.join(REPO_ROOT, 'tutorials', 'scripts')

DEFAULT_DATA_DIR = (
    '/Users/colby/Research/Programing/'
    'P3D_example_simulation_data/staging'
)
DATA_DIR = os.environ.get('PY3D_TUTORIAL_DATA', DEFAULT_DATA_DIR)


def _scripts():
    """Return absolute paths to the four Phase 6a tutorial scripts."""
    names = [
        '01_getting_started.py',
        '02_exploring_fields.py',
        '03_time_evolution.py',
        '04_pressure_tensor.py',
    ]
    return [os.path.join(TUTORIAL_SCRIPTS, n) for n in names]


pytestmark = pytest.mark.skipif(
    not os.path.isdir(DATA_DIR),
    reason=(f'Tutorial data dir not found: {DATA_DIR}. '
            f'Set PY3D_TUTORIAL_DATA to enable these tests.'),
)


@pytest.mark.slow
@pytest.mark.parametrize('script_path', _scripts(),
                         ids=lambda p: os.path.basename(p))
def test_tutorial_script(script_path, tmp_path):
    """Run one tutorial script end-to-end; require exit code 0."""
    env = os.environ.copy()
    env['MPLBACKEND'] = 'Agg'           # headless plotting
    env['PY3D_TUTORIAL_DATA'] = DATA_DIR
    # Make py3d importable in the subprocess. When pytest runs as `python
    # -m pytest`, py3d is on sys.path because cwd=repo root, but a
    # subprocess gets a fresh interpreter whose sys.path[0] is the
    # script's directory. Add the repo root explicitly.
    existing = env.get('PYTHONPATH', '')
    env['PYTHONPATH'] = (
        REPO_ROOT + (os.pathsep + existing if existing else '')
    )

    result = subprocess.run(
        [sys.executable, script_path],
        env=env,
        capture_output=True,
        text=True,
        cwd=REPO_ROOT,
        timeout=600,
    )

    assert result.returncode == 0, (
        f'{os.path.basename(script_path)} exited {result.returncode}\n'
        f'--- stdout ---\n{result.stdout}\n'
        f'--- stderr ---\n{result.stderr}'
    )
    # Every script prints at least one "Saved <png>" line.
    assert 'Saved' in result.stdout, (
        f'No "Saved" line in stdout — script may have skipped its plot.\n'
        f'stdout:\n{result.stdout}'
    )
