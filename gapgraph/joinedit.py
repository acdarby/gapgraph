"""Join Editor support (after Gap5's Join Editor).

Two contigs are laid out in one *join coordinate system*: contig A oriented
so the joining end is on its right and occupies [-len(A), 0); contig B
oriented so its joining end is on its left and starts at x = gap (negative
gap = overlap).  Reads and consensus from both are transformed into that
system so the two consensus lines sit next to each other in register,
separated by a differences line, and scroll together.
"""
from .seqs import revcomp

_COMP = str.maketrans("ACGTNacgtn*", "TGCANtgcan*")


def orient_read(r, L, reverse):
    """Transform a render_read() dict onto the reverse strand of a contig of length L."""
    if not reverse:
        return r
    s = r["s"].translate(_COMP)[::-1]
    x = L - (r["x"] + len(r["s"]))
    ins = [[L - p, revcomp(b) if b.isalpha() else b] for p, b in r["ins"]]
    out = dict(r)
    out.update(x=x, s=s, ins=ins, rev=not r["rev"])
    return out


def layout(cand_a, cand_b, gap, la, lb):
    """Return (A_reverse, B_reverse, A_offset, B_offset): join_x = offset + oriented_col."""
    a_rev = cand_a[1] == "L"      # A's joining end must be on the right
    b_rev = cand_b[1] == "R"      # B's joining end must be on the left
    return a_rev, b_rev, -la, gap


def banded_align(a, b, band=None, match=2, mismatch=-3, gap_open=-4):
    """Global alignment of a vs b within a diagonal band (linear gap cost).
    Returns (aligned_a, aligned_b) with '*' pads, Gap5 style."""
    n, m = len(a), len(b)
    if band is None:
        band = max(20, abs(n - m) + 20 + int(0.01 * max(n, m)))
    NEG = -10 ** 9
    W = 2 * band + 1
    # score/traceback rows indexed by j - i + band
    prev = [NEG] * W
    tb = []
    for d in range(W):
        j = d - band
        if 0 <= j <= m:
            prev[d] = gap_open * j
    tb.append([2] * W)  # row 0: all from left
    for i in range(1, n + 1):
        cur = [NEG] * W
        t = [0] * W
        ai = a[i - 1]
        for d in range(W):
            j = i + d - band
            if j < 0 or j > m:
                continue
            best, how = NEG, 0
            if j == 0:
                best, how = gap_open * i, 1
            else:
                # diagonal: (i-1, j-1) -> same d
                v = prev[d]
                if v > NEG:
                    v += match if ai == b[j - 1] else mismatch
                    best, how = v, 0
                # up: (i-1, j) -> d+1
                if d + 1 < W and prev[d + 1] > NEG and prev[d + 1] + gap_open > best:
                    best, how = prev[d + 1] + gap_open, 1
                # left: (i, j-1) -> d-1
                if d - 1 >= 0 and cur[d - 1] > NEG and cur[d - 1] + gap_open > best:
                    best, how = cur[d - 1] + gap_open, 2
            cur[d] = best
            t[d] = how
        tb.append(t)
        prev = cur
    # traceback from (n, m)
    i, j = n, m
    out_a, out_b = [], []
    while i > 0 or j > 0:
        d = j - i + band
        how = tb[i][d] if 0 <= d < W else (1 if i > 0 else 2)
        if i == 0:
            how = 2
        if j == 0:
            how = 1
        if how == 0:
            out_a.append(a[i - 1]); out_b.append(b[j - 1]); i -= 1; j -= 1
        elif how == 1:
            out_a.append(a[i - 1]); out_b.append("*"); i -= 1
        else:
            out_a.append("*"); out_b.append(b[j - 1]); j -= 1
    return "".join(reversed(out_a)), "".join(reversed(out_b))


def align_overlap(seq_a, seq_b, gap, max_len=5000):
    """Align the overlapping ends (A oriented end-right, B end-left).
    Returns a dict with identity and a Gap5-style text alignment."""
    if gap is None or gap >= 0:
        return {"overlap": 0, "gap": gap, "text": f"No overlap: {gap} bp gap between the ends "
                                                  f"(will be filled with Ns)."}
    ovl = -gap
    shown = min(ovl, max_len)
    # align the part of the overlap closest to the join point boundary
    a_part = seq_a[len(seq_a) - ovl: len(seq_a) - ovl + shown]
    b_part = seq_b[:shown]
    aa, bb = banded_align(a_part, b_part)
    cols = len(aa)
    matches = sum(1 for x, y in zip(aa, bb) if x == y and x != "*")
    mism = sum(1 for x, y in zip(aa, bb) if x != y and "*" not in (x, y))
    indel = sum(1 for x, y in zip(aa, bb) if "*" in (x, y))
    ident = 100.0 * matches / max(1, cols)
    lines = []
    for k in range(0, cols, 60):
        sa, sb = aa[k:k + 60], bb[k:k + 60]
        mid = "".join(":" if x == y else ("!" if "*" not in (x, y) else " ") for x, y in zip(sa, sb))
        lines += [f"A {sa}", f"  {mid}", f"B {sb}", ""]
    return {"overlap": ovl, "aligned": shown, "identity": round(ident, 2), "mismatches": mism,
            "indels": indel, "pct_mismatch": round(100 - ident, 2),
            "text": "\n".join(lines[:400]) + ("\n…(truncated)" if len(lines) > 400 else "")}
