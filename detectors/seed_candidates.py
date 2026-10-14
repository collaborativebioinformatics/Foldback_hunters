"""
numba candidate finder for `--mode seed` (requires numba + numpy).

Finds the first read position i whose k-mer matches a stride-k k-mer of
RC(read), plus the smallest q holding that k-mer. Same result as the pure-Python
dict scan in foldback_score.find_candidate_py, ~2.8x faster, and nogil so
`--threads` scales. Reads with characters outside ACGTN fall back to the
pure-Python scan so string-equality semantics stay identical.

Importing this module compiles the kernel (or loads the on-disk cache), so each
`--processes` worker pays that cost once at import.
"""

import numpy as np
from numba import njit

_LUT = np.full(256, 255, dtype=np.uint8)
for _c, _v in zip(b"ACGTN", (0, 1, 2, 3, 4)):
    _LUT[_c] = _v


def encode(seq):
    try:
        raw = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
    except UnicodeEncodeError:
        return None
    codes = _LUT[raw]
    if codes.max() == 255:
        return None
    return codes


@njit(nogil=True, cache=True)
def find_candidate_codes(codes, k):
    """Return (i, q) for encoded read `codes`; (-1, -1) if no k-mer matches."""
    n = codes.shape[0]
    nq = (n - k) // k + 1  # number of stride-k RC k-mers
    # open-addressing hash of base-5 k-mer keys -> smallest q
    size = 1
    while size < 2 * nq:
        size <<= 1
    mask = size - 1
    tkeys = np.full(size, -1, dtype=np.int64)
    tq = np.zeros(size, dtype=np.int64)
    for qi in range(nq):
        q = qi * k
        key = np.int64(0)
        for j in range(k):
            # RC code at position t = comp(codes[n-1-t]); comp(c) = 3-c for c<4 else 4
            c = np.int64(codes[n - 1 - (q + j)])
            if c < 4:
                c = 3 - c
            key = key * 5 + c
        h = ((key * 2654435761) >> 8) & mask
        while True:
            if tkeys[h] == -1:
                tkeys[h] = key
                tq[h] = q
                break
            if tkeys[h] == key:  # keep first (min) q
                break
            h = (h + 1) & mask
    top = np.int64(1)
    for _ in range(k - 1):
        top *= 5
    w = np.int64(0)
    for j in range(k):
        w = w * 5 + codes[j]
    i = 0
    last = n - k
    while True:
        h = ((w * 2654435761) >> 8) & mask
        while tkeys[h] != -1:
            if tkeys[h] == w:
                return i, tq[h]
            h = (h + 1) & mask
        if i == last:
            break
        w = (w - codes[i] * top) * 5 + codes[i + k]
        i += 1
    return -1, -1


def warmup():
    """Trigger cache load."""
    find_candidate_codes(np.zeros(64, dtype=np.uint8), 13)


warmup()
