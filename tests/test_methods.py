import numpy as np
import pytest

from py3d._methods import _convert, _num_to_ext, interp_field, load_param, vprint

# ---------------------------------------------------------------------------
# _convert
# ---------------------------------------------------------------------------

def test_convert_int():
    assert _convert("42") == 42
    assert isinstance(_convert("42"), int)


def test_convert_float():
    assert _convert("3.14") == pytest.approx(3.14)
    assert isinstance(_convert("3.14"), float)


def test_convert_string():
    result = _convert("hello")
    assert result == "hello"
    assert isinstance(result, str)


# ---------------------------------------------------------------------------
# _num_to_ext
# ---------------------------------------------------------------------------

def test_num_to_ext_zero():
    assert _num_to_ext(0) == "000"


def test_num_to_ext_single_digit():
    assert _num_to_ext(7) == "007"


def test_num_to_ext_triple_digit():
    assert _num_to_ext(999) == "999"


def test_num_to_ext_none():
    assert _num_to_ext(None) is None


# ---------------------------------------------------------------------------
# interp_field
# ---------------------------------------------------------------------------

def test_interp_field_2d_uniform():
    """Interpolating a constant 2D field should return the constant value."""
    fld = np.full((10, 10), 5.0)
    # r0 at the exact center of cell [5,5]: cell center = (i+0.5)*dl
    # dl = sim_lens/nl = [1.0, 1.0], cell 5 center = 5.5
    result = interp_field(fld, r0=[5.5, 5.5], sim_lens=[10., 10.])
    assert result == pytest.approx(5.0)


def test_interp_field_2d_midpoint():
    """At an exact cell center, result equals that cell's value."""
    # fld[i,j] = i + j; at cell [5,5], exact center r0=[5.5,5.5], wl=[0,0]
    fld = np.fromfunction(lambda i, j: i + j, (10, 10), dtype=float)
    result = interp_field(fld, r0=[5.5, 5.5], sim_lens=[10., 10.])
    assert result == pytest.approx(fld[5, 5])


def test_interp_field_3d_uniform():
    """Interpolating a constant 3D field should return the constant value."""
    fld = np.full((8, 8, 8), 3.0)
    result = interp_field(fld, r0=[4.5, 4.5, 4.5], sim_lens=[8., 8., 8.])
    assert result == pytest.approx(3.0)


# ---------------------------------------------------------------------------
# load_param
# ---------------------------------------------------------------------------

def test_load_param_parses_ints(param_file):
    _, fpath = param_file
    param = load_param(param_file=fpath)
    assert param['nx'] == 64
    assert param['ny'] == 32
    assert isinstance(param['nx'], int)


def test_load_param_parses_floats(param_file):
    _, fpath = param_file
    param = load_param(param_file=fpath)
    assert param['lx'] == pytest.approx(10.0)
    assert isinstance(param['lx'], float)


def test_load_param_stores_filename(param_file):
    _, fpath = param_file
    param = load_param(param_file=fpath)
    assert 'file' in param


# ---------------------------------------------------------------------------
# vprint
# ---------------------------------------------------------------------------

class _Verbose:
    _verbose = True


class _Silent:
    _verbose = False


def test_vprint_prints_when_verbose(capsys):
    vprint(_Verbose(), "hello", " world")
    assert "hello world" in capsys.readouterr().out


def test_vprint_silent(capsys):
    vprint(_Silent(), "should not appear")
    assert capsys.readouterr().out == ""


def test_vprint_no_attribute(capsys):
    """Objects without _verbose should behave as silent."""
    vprint(object(), "should not appear")
    assert capsys.readouterr().out == ""


# ---------------------------------------------------------------------------
# Deferred: Movie, Dump, DumpID
# ---------------------------------------------------------------------------
# These classes use input() prompts and require binary simulation files.
# They can be tested once Phase 5 refactors them to accept explicit paths.

@pytest.mark.skip(reason="Movie requires binary sim files and input() prompts; revisit in Phase 5")
def test_movie_init_stub():
    pass


@pytest.mark.skip(reason="Dump requires binary sim files and input() prompts; revisit in Phase 5")
def test_dump_init_stub():
    pass
