import numpy as np
import pytest


@pytest.fixture
def simple_2d_field():
    """Minimal field dict: 20x10 grid with bx=sin(x), by=cos(y)."""
    xx = np.linspace(0, 10, 20)
    yy = np.linspace(0, 5, 10)
    bx = np.sin(xx[:, None]) * np.ones((20, 10))
    by = np.ones((20, 10)) * np.cos(yy)
    return {'xx': xx, 'yy': yy, 'bx': bx, 'by': by}


@pytest.fixture
def param_file(tmp_path):
    """Minimal P3D .in param file written to a temp directory."""
    content = "\n".join([
        "#define nprocs 4",
        "#define nx 64",
        "#define ny 32",
        "#define nz 1",
        "#define lx 10.0",
        "#define ly 5.0",
        "#define prk 4",
    ])
    f = tmp_path / "p3d.in"
    f.write_text(content)
    return str(tmp_path), str(f)


@pytest.fixture
def pressure_tensor_field():
    """Field dict with all components needed by rotate_ten (default path)."""
    rng = np.random.default_rng(42)
    shape = (20, 10)
    d = {}
    # B-field (unit z-hat so parallel = pizzav)
    d['bxav'] = np.zeros(shape)
    d['byav'] = np.zeros(shape)
    d['bzav'] = np.ones(shape)
    # Diagonal pressure tensor
    d['pixxav'] = rng.random(shape)
    d['piyyav'] = rng.random(shape)
    d['pizzav'] = rng.random(shape)
    # Off-diagonal (zero for simplicity)
    d['pixyav'] = np.zeros(shape)
    d['piyzav'] = np.zeros(shape)
    d['pixzav'] = np.zeros(shape)
    return d
