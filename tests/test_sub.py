import numpy as np
import pytest

from py3d.sub import calc_pdf, calc_psi, findval, rotate_ten, var_at

# ---------------------------------------------------------------------------
# findval
# ---------------------------------------------------------------------------

def test_findval_exact_match():
    vec = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    assert findval(vec, 3.0) == 2


def test_findval_nearest_below():
    vec = np.linspace(0.0, 10.0, 11)   # [0, 1, 2, ..., 10]
    assert findval(vec, 2.4) == 2       # closest to index 2 (value 2.0)


def test_findval_nearest_above():
    vec = np.linspace(0.0, 10.0, 11)
    assert findval(vec, 2.6) == 3       # closest to index 3 (value 3.0)


def test_findval_first_element():
    vec = np.linspace(0.0, 10.0, 11)
    assert findval(vec, 0.0) == 0


def test_findval_last_element():
    vec = np.linspace(0.0, 10.0, 11)
    assert findval(vec, 10.0) == 10


# ---------------------------------------------------------------------------
# calc_pdf
# ---------------------------------------------------------------------------

def test_calc_pdf_returns_two_arrays():
    ar = np.random.default_rng(0).standard_normal(1000)
    result = calc_pdf(ar)
    assert result is not None
    binvals, pdf = result
    assert len(binvals) == len(pdf)


def test_calc_pdf_output_length():
    """Number of bins = len(ar) // weight."""
    ar = np.ones(500)
    ar = np.random.default_rng(1).standard_normal(500)
    binvals, pdf = calc_pdf(ar, weight=50)
    assert len(binvals) == 500 // 50


def test_calc_pdf_empty_array_returns_none():
    result = calc_pdf(np.array([]))
    assert result is None


def test_calc_pdf_inc_path_runs():
    """inc > 0 takes the increment branch without error."""
    ar = np.random.default_rng(2).standard_normal((10, 10))
    result = calc_pdf(ar, inc=1)
    assert result is not None


# ---------------------------------------------------------------------------
# calc_psi
# ---------------------------------------------------------------------------

def test_calc_psi_shape(simple_2d_field):
    psi = calc_psi(simple_2d_field)
    assert psi.shape == simple_2d_field['bx'].shape


def test_calc_psi_uniform_bx():
    """Uniform bx=1, by=0 → psi linear in y, constant in x."""
    nx, ny = 20, 10
    xx = np.linspace(0, 10, nx)
    yy = np.linspace(0, 5, ny)
    dy = yy[1] - yy[0]
    d = {
        'xx': xx, 'yy': yy,
        'bx': np.ones((nx, ny)),
        'by': np.zeros((nx, ny)),
    }
    psi = calc_psi(d)
    expected_row = np.arange(ny) * dy
    assert np.allclose(psi[0, :], expected_row)
    assert np.allclose(psi[5, :], psi[0, :])


def test_calc_psi_uses_bxav_when_present():
    """Uses bxav/byav keys when available."""
    nx, ny = 10, 8
    xx = np.linspace(0, 5, nx)
    yy = np.linspace(0, 4, ny)
    d = {
        'xx': xx, 'yy': yy,
        'bx': np.zeros((nx, ny)),   # would give psi=0
        'by': np.zeros((nx, ny)),
        'bxav': np.ones((nx, ny)),  # uniform bxav → non-zero psi
        'byav': np.zeros((nx, ny)),
    }
    psi = calc_psi(d)
    assert not np.all(psi == 0)


# ---------------------------------------------------------------------------
# var_at
# ---------------------------------------------------------------------------

def test_var_at_constant_field():
    """Interpolating a constant field should return that constant."""
    xx = np.arange(10, dtype=float)
    yy = np.arange(10, dtype=float)
    fdic = {
        'xx': xx,
        'yy': yy,
        'f': np.full((10, 10), 7.0),
    }
    result = var_at(fdic, 'f', r0=[3.0, 4.0])
    assert result == pytest.approx(7.0)


def test_var_at_result_in_range():
    """Interpolated value should lie between min and max of field."""
    rng = np.random.default_rng(5)
    xx = np.linspace(0, 9, 20)
    yy = np.linspace(0, 9, 20)
    field = rng.random((20, 20))
    fdic = {'xx': xx, 'yy': yy, 'f': field}
    result = var_at(fdic, 'f', r0=[4.5, 4.5])
    assert field.min() <= result <= field.max()


# ---------------------------------------------------------------------------
# rotate_ten
# ---------------------------------------------------------------------------

def test_rotate_ten_adds_par_perp_keys(pressure_tensor_field):
    d = pressure_tensor_field.copy()
    rotate_ten(d)
    assert 'piparav' in d
    assert 'piperp1av' in d
    assert 'piperp2av' in d


def test_rotate_ten_output_shape(pressure_tensor_field):
    d = pressure_tensor_field.copy()
    rotate_ten(d)
    expected_shape = pressure_tensor_field['pixxav'].shape
    assert d['piparav'].shape == expected_shape


def test_rotate_ten_z_aligned_parallel_equals_pizz(pressure_tensor_field):
    """With B along z-hat, pipar should equal pizz."""
    d = pressure_tensor_field.copy()
    rotate_ten(d)
    assert np.allclose(d['piparav'], d['pizzav'])


def test_rotate_ten_skips_if_already_present(pressure_tensor_field):
    """If output key already exists, rotate_ten should leave it unchanged."""
    d = pressure_tensor_field.copy()
    sentinel = np.full_like(d['pixxav'], 999.0)
    d['piparav'] = sentinel.copy()
    rotate_ten(d, overwrite=False)
    assert np.all(d['piparav'] == 999.0)


# ---------------------------------------------------------------------------
# Deferred: plotting functions
# ---------------------------------------------------------------------------
# ims, PatPlotter.make_plots(), VDistPlotter.plot2d() produce matplotlib
# figures. Automated testing requires visual regression tooling (e.g.
# pytest-mpl) which is out of scope for this phase.

@pytest.mark.skip(reason="Plotting tests require visual regression tooling; revisit post-Phase 4")
def test_ims_stub():
    pass
