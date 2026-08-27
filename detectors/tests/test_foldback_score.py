"""
Unit tests for detectors/foldback_score.py.

Fixtures build foldback reads (read = A + RC(B), fold point p = len(A),
tail t = len(B)) with construction math reimplemented locally, independent
of sim/simulate_foldbacks_v2.py.
"""

import random

import pytest
from foldback_score import ScoreResult, revcomp, score_read_full, score_read_probe

DEFAULT_PROBE_LENGTHS = (50, 150, 500, 1500)


def make_foldback_read(L, p, seed=0):
    """Build a noiseless foldback read: A + RC(B), true fold point p."""
    rng = random.Random(seed)
    seq = "".join(rng.choice("ACGT") for _ in range(L))
    t = min(p, L - p)
    first = seq[:p]
    tail_source = first[-t:]
    read = first + revcomp(tail_source)
    assert len(read) == p + t
    return read, t


# ---------------------------------------------------------------------------
# probe mode: position inference
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("L,p", [
    (5000, 2500),   # middle: t=2500, well above the 1500 probe-length ceiling
    (5000, 3750),   # off_center: t=1250
    (5000, 4850),   # near_end: t=150, so only k=50 and k=150 qualify (<= 2t=300)
])
def test_probe_position_exact(L, p):
    read, _t = make_foldback_read(L, p, seed=1)
    result = score_read_probe("r1", read, probe_lengths=DEFAULT_PROBE_LENGTHS, min_len=1000)
    assert result.status == "ok"
    assert result.raw_score == 1.0
    assert result.fold_position_bp == p


# ---------------------------------------------------------------------------
# full mode: position inference (exact regardless of t)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("L,p", [
    (5000, 2500),
    (5000, 3750),
    pytest.param(5000, 4850, marks=pytest.mark.xfail(
        reason="Real finding: for small t (near_end reads), a longer/gappier "
               "chance alignment in the random background can out-score the "
               "true fold signal in full mode's max-score search.",
        strict=True,
    )),
])
def test_full_position_exact(L, p):
    read, _t = make_foldback_read(L, p, seed=5)
    result = score_read_full("r5", read, min_len=1000)
    assert result.status == "ok"
    assert result.raw_score == pytest.approx(1.0)
    assert result.fold_position_bp == p


# ---------------------------------------------------------------------------
# too-short guard (both modes, never raises)
# ---------------------------------------------------------------------------

def test_too_short_probe():
    seq = "ACGT" * 125  # 500bp < default min_len 1000
    result = score_read_probe("short1", seq)
    assert result == ScoreResult("short1", None, None, "too_short")


def test_too_short_full():
    seq = "ACGT" * 125  # 500bp < default min_len 1000
    result = score_read_full("short2", seq)
    assert result == ScoreResult("short2", None, None, "too_short")


# ---------------------------------------------------------------------------
# N-base handling (no special-casing; N is an ordinary mismatching symbol)
# ---------------------------------------------------------------------------

def test_n_base_full_pinned_score():
    # p = t = 100 -> whole 200bp read is one ungapped perfect-match span.
    # A single mutated base appears on both sides of read-vs-RC(read), so
    # identity drops by 2 bases, not 1.
    L, p = 200, 100
    read, t = make_foldback_read(L, p, seed=6)
    assert t == 100
    mutated = read[:100] + "N" + read[101:]
    result = score_read_full("n1", mutated, min_len=10)
    assert result.status == "ok"
    assert result.raw_score == pytest.approx(198 / 200)
    assert abs(result.fold_position_bp - p) <= 1


def test_n_base_probe_no_special_casing():
    # N placed inside the probe's matched window; edlib treats it as an
    # ordinary mismatch (editDistance=1, score=(k-1)/k).
    L, p = 200, 100
    read, t = make_foldback_read(L, p, seed=7)
    assert t == 100
    mutated = read[:10] + "N" + read[11:]
    result = score_read_probe("n2", mutated, probe_lengths=(50,), min_len=10)
    assert result.status == "ok"
    assert result.raw_score == pytest.approx(49 / 50)
