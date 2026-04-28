import numpy as np
import pytest

from py3d.vdist import VDist


@pytest.fixture
def vdist():
    return VDist()


@pytest.fixture
def particle_velocities():
    rng = np.random.default_rng(42)
    n = 2000
    return {
        'v1': rng.standard_normal(n),
        'v2': rng.standard_normal(n),
        'v3': rng.standard_normal(n),
    }


# ---------------------------------------------------------------------------
# vdist2d
# ---------------------------------------------------------------------------

def test_vdist2d_output_shapes(vdist, particle_velocities):
    v = particle_velocities
    H, xe, ye = vdist.vdist2d(v['v1'], v['v2'], bins=20, range=[[-5, 5], [-5, 5]])
    assert H.shape == (20, 20)
    assert xe.shape == (21,)
    assert ye.shape == (21,)


def test_vdist2d_all_particles_counted(vdist):
    """With a range that covers all particles, H sums to N."""
    rng = np.random.default_rng(7)
    v1 = rng.uniform(-1, 1, 500)
    v2 = rng.uniform(-1, 1, 500)
    H, _, _ = vdist.vdist2d(v1, v2, bins=10, range=[[-1, 1], [-1, 1]])
    assert H.sum() == 500


def test_vdist2d_uniform_weights(vdist):
    """Uniform weight=w scales histogram by w vs. weight=1."""
    rng = np.random.default_rng(9)
    v1 = rng.standard_normal(200)
    v2 = rng.standard_normal(200)
    kwargs = dict(bins=10, range=[[-4, 4], [-4, 4]])
    H_plain, _, _ = vdist.vdist2d(v1, v2, **kwargs)
    H_weighted, _, _ = vdist.vdist2d(v1, v2, weights=np.full(200, 2.0), **kwargs)
    assert np.allclose(H_weighted, H_plain * 2.0)


def test_vdist2d_v3_trim_reduces_count(vdist, particle_velocities):
    """Providing v3 with a tight dz should drop particles."""
    v = particle_velocities
    # dz much smaller than std of v3 → many particles trimmed
    H_trim, _, _ = vdist.vdist2d(v['v1'], v['v2'], v3=v['v3'], dz=0.1,
                                  bins=10, range=[[-4, 4], [-4, 4]])
    H_full, _, _ = vdist.vdist2d(v['v1'], v['v2'],
                                  bins=10, range=[[-4, 4], [-4, 4]])
    assert H_trim.sum() < H_full.sum()


# ---------------------------------------------------------------------------
# vdist2d_pitch
# ---------------------------------------------------------------------------

def test_vdist2d_pitch_output_shapes(vdist, particle_velocities):
    v = particle_velocities
    H, xe, ye = vdist.vdist2d_pitch(
        v['v1'], v['v2'], v['v3'], pa=90., dpa=20.,
        bins=15, range=[[-5, 5], [-5, 5]]
    )
    assert H.shape == (15, 15)
    assert xe.shape == (16,)
    assert ye.shape == (16,)


def test_vdist2d_pitch_tight_filter_returns_fewer(vdist, particle_velocities):
    """A tight pitch-angle window selects fewer particles than a wide one."""
    v = particle_velocities
    kwargs = dict(bins=10, range=[[-5, 5], [-5, 5]])
    H_wide, _, _ = vdist.vdist2d_pitch(v['v1'], v['v2'], v['v3'],
                                        pa=90., dpa=180., **kwargs)
    H_tight, _, _ = vdist.vdist2d_pitch(v['v1'], v['v2'], v['v3'],
                                         pa=90., dpa=5., **kwargs)
    # Pre-normalisation counts: count non-zero bins as proxy
    assert np.count_nonzero(H_tight) <= np.count_nonzero(H_wide)


def test_vdist2d_pitch_different_angles_differ(vdist, particle_velocities):
    """Selecting pa=0 vs pa=90 on isotropic data should give different maps."""
    v = particle_velocities
    kwargs = dict(bins=10, range=[[-5, 5], [-5, 5]])
    H_para, _, _ = vdist.vdist2d_pitch(v['v1'], v['v2'], v['v3'],
                                        pa=0., dpa=20., **kwargs)
    H_perp, _, _ = vdist.vdist2d_pitch(v['v1'], v['v2'], v['v3'],
                                        pa=90., dpa=20., **kwargs)
    assert not np.allclose(H_para, H_perp)


# ---------------------------------------------------------------------------
# spec1d
# ---------------------------------------------------------------------------

def test_spec1d_shape(vdist, particle_velocities):
    v = particle_velocities
    n = len(v['v1'])
    pts = {
        'x': np.random.default_rng(0).uniform(0, 10, n),
        'y': np.random.default_rng(1).uniform(0, 10, n),
        'z': np.random.default_rng(2).uniform(0, 10, n),
        'v0': v['v1'],
        'v1': v['v2'],
        'v2': v['v3'],
    }
    H, xe, ye = vdist.spec1d(pts, dir='x', pa=90., dpa=45., mass=1.0, bins=10)
    assert xe.shape[0] == ye.shape[0]


def test_spec1d_all_directions(vdist, particle_velocities):
    """spec1d should run for dir='x', 'y', 'z' without error."""
    v = particle_velocities
    n = len(v['v1'])
    rng = np.random.default_rng(3)
    pts = {
        'x': rng.uniform(0, 10, n), 'y': rng.uniform(0, 10, n),
        'z': rng.uniform(0, 10, n),
        'v0': v['v1'], 'v1': v['v2'], 'v2': v['v3'],
    }
    for d in ('x', 'y', 'z'):
        vdist.spec1d(pts, dir=d, pa=90., dpa=45., mass=1.0, bins=10)


# ---------------------------------------------------------------------------
# Deferred: PartTrace TPRun
# ---------------------------------------------------------------------------
# TPRun requires compiled C extensions (functions.so / functions_rel.so).
# Tests should be added after confirming a build step in CI (Phase 4).

@pytest.mark.skip(reason="TPRun tests require compiled C extensions; add build step in Phase 4 CI")
def test_tprun_stub():
    pass
