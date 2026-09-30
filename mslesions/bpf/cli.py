"""bpf T1w brainmask [lesions] | bpf --manifest m.tsv (subject, t1, brainmask[, mask][, age])"""

import argparse
import json
import sys

import pandas as pd

from .core import bpf


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="bpf", description=__doc__)
    ap.add_argument("t1", nargs="?")
    ap.add_argument("brain", nargs="?")
    ap.add_argument("lesions", nargs="?")
    ap.add_argument("--manifest")
    a = ap.parse_args(argv)
    if a.manifest:
        m = pd.read_csv(a.manifest, sep="\t")
        from msdataqc.tissue import TissueError

        rows = []
        for _, r in m.iterrows():
            try:
                rows.append({**r.to_dict(), **bpf(r["t1"], r["brainmask"], r.get("mask"))})
            except TissueError as e:
                print(r["subject"], "SKIPPED:", e, file=sys.stderr)
        out = pd.DataFrame(rows)
        out.drop(columns=[c for c in ("t1", "brainmask", "mask") if c in out]).to_csv(sys.stdout, sep="\t", index=False)
        if "age" in out:
            print(f"Spearman BPF vs age: {out['bpf'].corr(out['age'], method='spearman'):.2f} (n={len(out)})", file=sys.stderr)
        return 0
    print(json.dumps(bpf(a.t1, a.brain, a.lesions), indent=2))
    return 0
