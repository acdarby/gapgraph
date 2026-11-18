"""Sequence helpers and FASTA I/O."""
import gzip

_COMP = str.maketrans("ACGTNacgtnRYKMSWBDHVrykmswbdhv", "TGCANtgcanYRMKSWVHDByrmkswvhdb")


def revcomp(s):
    return s.translate(_COMP)[::-1]


def open_text(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line)
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def write_fasta(records, path, width=80):
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


def kmer_index(seq, k):
    """Map k-mer -> list of positions (uppercase, N-containing k-mers skipped)."""
    idx = {}
    s = seq.upper()
    for i in range(len(s) - k + 1):
        km = s[i:i + k]
        if "N" in km:
            continue
        idx.setdefault(km, []).append(i)
    return idx
