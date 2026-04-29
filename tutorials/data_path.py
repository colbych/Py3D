#######################################################################
#                                                                     #
#                  Python Progs :  data_path.py                       #
#                  Aruthor      :  Colby Haggerty                     #
#                  Date         :  2026.04.28                         #
#                                                                     #
#######################################################################
"""
Single source of truth for the tutorial dataset path.

All tutorial scripts and notebooks import DATA_DIR and PARAM_FILE from
here so a user with the example simulation in a different location can
override it once via the PY3D_TUTORIAL_DATA environment variable:

    export PY3D_TUTORIAL_DATA=/path/to/your/staging
"""

import os


DEFAULT_DATA_DIR = (
    '/Users/colby/Research/Programing/'
    'P3D_example_simulation_data/staging'
)

DATA_DIR = os.environ.get('PY3D_TUTORIAL_DATA', DEFAULT_DATA_DIR)
PARAM_FILE = 'param_asym00p'
NAME_STYLE = 'p3d'   # files are 'movie.bx.000', not 'bx'


def require_data_dir():
    """Raise FileNotFoundError if DATA_DIR is missing.

    Tutorials call this at the top so failure is immediate and points the
    user at the override env var instead of failing inside py3d code.
    """
    if not os.path.isdir(DATA_DIR):
        raise FileNotFoundError(
            f'Tutorial data directory not found: {DATA_DIR}\n'
            f'Set PY3D_TUTORIAL_DATA to point at your staging/ folder, e.g.:\n'
            f'    export PY3D_TUTORIAL_DATA=/path/to/staging'
        )
