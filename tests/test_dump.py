import pytest

from py3d.dump import Dump
from py3d.dumpID import DumpID


# ---------------------------------------------------------------------------
# Error paths — interactive=False must raise instead of prompting on stdin.
# ---------------------------------------------------------------------------

def test_dump_no_files_raises(tmp_path, param_file):
    """Empty directory + interactive=False -> FileNotFoundError, no prompt."""
    _, fpath = param_file
    with pytest.raises(FileNotFoundError):
        Dump(num=0, param_file=fpath, path=str(tmp_path), interactive=False)


def test_dump_invalid_num_raises(tmp_path, param_file):
    """An existing dump-000 with num=99 -> ValueError listing valid choices."""
    _, fpath = param_file
    (tmp_path / "p3d-001.000").touch()
    with pytest.raises(ValueError, match="Valid choices"):
        Dump(num=99, param_file=fpath, path=str(tmp_path), interactive=False)


def test_dumpid_no_files_raises(tmp_path, param_file):
    _, fpath = param_file
    with pytest.raises(FileNotFoundError):
        DumpID(num=0, param_file=fpath, path=str(tmp_path), interactive=False)


# ---------------------------------------------------------------------------
# Happy path — minimal scaffolding (empty p3d-001.000) is enough for the
# constructor; read_particles/read_fields require real binary content and
# are validated against real simulation data outside the repo (Phase 5d).
# ---------------------------------------------------------------------------

def test_dump_constructs_with_explicit_args(tmp_path, param_file):
    _, fpath = param_file
    (tmp_path / "p3d-001.000").touch()
    d = Dump(num=0, param_file=fpath, path=str(tmp_path), interactive=False)
    assert d.num == "000"
    assert d.path == str(tmp_path)
    assert d.param["prk"] == 4


def test_dumpid_constructs_with_explicit_args(tmp_path, param_file):
    _, fpath = param_file
    (tmp_path / "p3d-001.000").touch()
    did = DumpID(num=0, param_file=fpath, path=str(tmp_path), interactive=False)
    assert did.dump.num == "000"
    assert did.param["prk"] == 4
