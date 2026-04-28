import pytest

from py3d.movie import Movie


# ---------------------------------------------------------------------------
# Error paths — interactive=False must raise instead of prompting on stdin.
# Happy-path construction needs real movie binary + log files and is
# validated against real simulation data outside the repo (Phase 5d).
# ---------------------------------------------------------------------------

def test_movie_no_files_raises(tmp_path, param_file):
    """Empty directory + interactive=False -> FileNotFoundError, no prompt."""
    _, fpath = param_file
    with pytest.raises(FileNotFoundError):
        Movie(num=0, param_file=fpath, path=str(tmp_path), interactive=False)


def test_movie_invalid_num_raises(tmp_path, param_file):
    """Two log files exist but num is bogus -> ValueError listing choices."""
    _, fpath = param_file
    (tmp_path / "movie.log.000").touch()
    (tmp_path / "movie.log.001").touch()
    with pytest.raises(ValueError, match="Valid choices"):
        Movie(num=99, param_file=fpath, path=str(tmp_path), interactive=False)
