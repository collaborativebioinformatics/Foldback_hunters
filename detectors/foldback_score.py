"""
Reference-free foldback read-level scoring (Task C).

Signature: read ~= forward_part + RC(segment_near_fold).

Two modes:
  - probe: edlib infix search of RC(tail) within the read.
  - full:  parasail Smith-Waterman of read vs RC(read).

revcomp() is reimplemented locally (does not import sim/simulate_foldbacks_v2.py).
"""

from dataclasses import dataclass
from typing import Optional

import edlib
import parasail

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")

DEFAULT_PROBE_LENGTHS = (50, 150, 500, 1500)
DEFAULT_MIN_LEN = 1000

# EMBOSS/dnafull-style gap penalties for full mode (confirmed with user).
GAP_OPEN = 10
GAP_EXTEND = 1


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


@dataclass
class ScoreResult:
    read_id: str
    raw_score: Optional[float]
    fold_position_bp: Optional[int]
    status: str


def score_read_probe(read_id, seq, probe_lengths=DEFAULT_PROBE_LENGTHS,
                      min_len=DEFAULT_MIN_LEN) -> ScoreResult:
    if len(seq) < min_len:
        return ScoreResult(read_id, None, None, "too_short")

    seq_rc = revcomp(seq)
    best_score = None
    best_k = None
    best_hit = None

    for k in probe_lengths:
        k = min(k, len(seq))
        probe = seq_rc[:k]
        result = edlib.align(probe, seq, mode="HW", task="locations")
        edit_distance = result["editDistance"]
        if edit_distance < 0:
            continue
        score = 1 - edit_distance / k
        # On ties, prefer the larger k (more specific match).
        if best_score is None or score > best_score or (
            score == best_score and k > best_k
        ):
            best_score = score
            best_k = k
            locations = result["locations"]
            best_hit = locations[0] if locations else None

    if best_score is None or best_hit is None:
        return ScoreResult(read_id, 0.0, None, "no_match")

    hit_start, hit_end = best_hit
    # For read = A + RC(B), fold point p = len(A), tail t = len(B):
    # hit_start == p - t for any k in (0, 2t], and len(seq) == p + t, so
    # averaging the two cancels t and recovers p exactly (k-independent).
    fold_position_bp = (hit_start + len(seq)) // 2

    return ScoreResult(read_id, best_score, fold_position_bp, "ok")


def score_read_full(read_id, seq, min_len=DEFAULT_MIN_LEN) -> ScoreResult:
    if len(seq) < min_len:
        return ScoreResult(read_id, None, None, "too_short")

    seq_rc = revcomp(seq)
    result = parasail.sw_trace(seq, seq_rc, GAP_OPEN, GAP_EXTEND, parasail.dnafull)

    comp = result.traceback.comp
    aligned_len = len(comp)
    if aligned_len == 0:
        return ScoreResult(read_id, 0.0, None, "no_match")

    matches = comp.count("|")
    identity = matches / aligned_len

    query_traceback = result.traceback.query
    matched_len_query = len(query_traceback) - query_traceback.count("-")
    start_query_0idx = result.end_query - matched_len_query + 1
    fold_position_bp = (start_query_0idx + result.end_query) // 2 + 1

    return ScoreResult(read_id, identity, fold_position_bp, "ok")


def score_read(read_id, seq, mode="probe", **kwargs) -> ScoreResult:
    if mode == "probe":
        return score_read_probe(read_id, seq, **kwargs)
    elif mode == "full":
        return score_read_full(read_id, seq, **kwargs)
    else:
        raise ValueError(f"Unknown mode: {mode!r} (expected 'probe' or 'full')")
