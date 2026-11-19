"""Aligning soft-clipped read ends ("hidden data" in Gap5 terms).

A soft clip is anchored at the column where its read stopped aligning, so the
question is always the same: how does this clip align to the sequence that
lies beyond the anchor?  `prefix_align` answers it with the query start fixed
at the target start and both ends free; the alignment is cut where it stops
being supported (X-drop), so adapter or junk tails stay clipped.

edlib (C++, MIT licence) does the alignment when installed; a banded pure
Python DP is the fallback.
"""
from collections import Counter

try:
    import edlib  # noqa
    HAVE_EDLIB = True
except Exception:  # pragma: no cover
    HAVE_EDLIB = False

import re

_EDCIG = re.compile(r"(\d+)([=XID])")


def _xdrop_trim(cols, match=1, mismatch=-2, gap=-2, xdrop=30):
    """cols: list of '=', 'X', 'I', 'D' per alignment column.
    Return the number of columns to keep (prefix with the best score)."""
    s = best = 0
    keep = 0
    for i, c in enumerate(cols):
        s += match if c == "=" else (mismatch if c == "X" else gap)
        if s > best:
            best, keep = s, i + 1
        elif best - s > xdrop:
            break
    return keep


def _cols_edlib(q, t):
    r = edlib.align(q, t, mode="SHW", task="path")
    cols = []
    for n, op in _EDCIG.findall(r["cigar"] or ""):
        cols.extend(op * int(n))
    return cols


def _cols_python(q, t, band=40, max_len=3000):
    q, t = q[:max_len], t[:max_len + band]
    n, m = len(q), len(t)
    NEG = -10 ** 9
    W = 2 * band + 1
    score = [[NEG] * W for _ in range(n + 1)]
    tb = [[0] * W for _ in range(n + 1)]
    for d in range(band, W):
        j = d - band
        if j <= m:
            score[0][d] = -2 * j
            tb[0][d] = 2
    best, bi, bd = 0, 0, band
    for i in range(1, n + 1):
        for d in range(W):
            j = i + d - band
            if j < 0 or j > m:
                continue
            v, how = NEG, 0
            if j > 0 and score[i - 1][d] > NEG:
                v = score[i - 1][d] + (1 if q[i - 1] == t[j - 1] else -2)
            if d + 1 < W and score[i - 1][d + 1] > NEG and score[i - 1][d + 1] - 2 > v:
                v, how = score[i - 1][d + 1] - 2, 1
            if d - 1 >= 0 and score[i][d - 1] > NEG and score[i][d - 1] - 2 > v:
                v, how = score[i][d - 1] - 2, 2
            if j == 0 and v == NEG:
                v, how = -2 * i, 1
            score[i][d], tb[i][d] = v, how
            if v > best:
                best, bi, bd = v, i, d
    cols = []
    i, d = bi, bd
    while i > 0 or d != band:
        j = i + d - band
        if i == 0:
            cols.append("D"); d -= 1; continue
        how = tb[i][d]
        if how == 0:
            cols.append("=" if q[i - 1] == t[j - 1] else "X"); i -= 1
        elif how == 1:
            cols.append("I"); i -= 1; d += 1
        else:
            cols.append("D"); d -= 1
    return cols[::-1]


def prefix_align(query, target):
    """Align `query` (anchored at its first base) to the start of `target`.

    Returns {"cols": [...], "q": query bases used, "t": target bases used,
    "identity": fraction of '=' columns} after X-drop trimming."""
    q, t = query.upper(), target.upper()
    if not q or not t:
        return {"cols": [], "q": 0, "t": 0, "identity": 0.0}
    cols = _cols_edlib(q, t) if HAVE_EDLIB else _cols_python(q, t)
    cols = cols[:_xdrop_trim(cols)]
    while cols and cols[-1] in "ID":
        cols.pop()
    nq = sum(1 for c in cols if c != "D")
    nt = sum(1 for c in cols if c != "I")
    ident = sum(1 for c in cols if c == "=") / max(1, len(cols))
    return {"cols": cols, "q": nq, "t": nt, "identity": ident}


def cols_to_cigar(cols):
    """Alignment columns -> SAM cigar ops [(n, op)] with = / X merged to M."""
    out = []
    for c in cols:
        op = "M" if c in "=X" else c
        if out and out[-1][1] == op:
            out[-1][0] += 1
        else:
            out.append([1, op])
    return [(n, op) for n, op in out]


def render_cols(query, cols):
    """Per-target-column characters ('*' = deletion) and insertions
    [(target_col, bases)] for an aligned query."""
    chars, ins = [], []
    qi = 0
    pending = ""
    for c in cols:
        if c == "I":
            pending += query[qi]
            qi += 1
            continue
        if pending:
            ins.append((len(chars), pending))
            pending = ""
        if c == "D":
            chars.append("*")
        else:
            chars.append(query[qi])
            qi += 1
    return "".join(chars), ins


def clip_consensus(clips, rounds=2, max_len=4000, min_depth=3):
    """Consensus of outward-oriented clips by iterated star alignment.
    Returns the consensus string (outward orientation)."""
    clips = [c.upper()[:max_len] for c in clips if len(c) >= 20]
    if len(clips) < min_depth:
        return ""
    lens = sorted(len(c) for c in clips)
    seed_len = lens[int(0.75 * (len(lens) - 1))]
    cons = next(c for c in sorted(clips, key=lambda c: abs(len(c) - seed_len)))
    for _ in range(rounds):
        cols_votes = [Counter() for _ in range(len(cons))]
        ins_votes = [Counter() for _ in range(len(cons) + 1)]
        cover = [0] * (len(cons) + 1)
        for c in clips:
            a = prefix_align(c, cons)
            chars, ins = render_cols(c[:a["q"]], a["cols"])
            for j, ch in enumerate(chars):
                cols_votes[j][ch] += 1
                cover[j] += 1
            for j, b in ins:
                ins_votes[j][b] += 1
        out = []
        for j in range(len(cons)):
            if cover[j] < min_depth:
                break
            iv = ins_votes[j]
            if iv and sum(iv.values()) * 2 > cover[j]:
                out.append(iv.most_common(1)[0][0])
            b = cols_votes[j].most_common(1)[0][0]
            if b != "*":
                out.append(b)
        new = "".join(out)
        if new == cons:
            break
        cons = new
    return cons
