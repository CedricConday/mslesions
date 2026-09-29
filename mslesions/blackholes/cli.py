"""blackholes T1w brainmask lesions [--lesions-out t.tsv] | blackholes --manifest m.tsv (subject, t1, brainmask, mask)"""

import argparse
import json
import sys

import pandas as pd

from .core import detect


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="blackholes", description=__doc__)
    ap.add_argument("t1", nargs="?")
    ap.add_argument("brain", nargs="?")
    ap.add_argument("lesions", nargs="?")
    ap.add_argument("--manifest")
    ap.add_argument("--lesions-out")
    a = ap.parse_args(argv)
    if a.manifest:
        rows = []
        from msdataqc.tissue import TissueError

        for _, r in pd.read_csv(a.manifest, sep="\t").iterrows():
            try:
                _, s = detect(r["t1"], r["brainmask"], r["mask"])
            except TissueError as e:
                rows.append({"subject": r["subject"], "error": str(e)})
                print(r["subject"], "SKIPPED:", e, file=sys.stderr, flush=True)
                continue
            rows.append({"subject": r["subject"], **s})
            print(r["subject"], file=sys.stderr, flush=True)
        pd.DataFrame(rows).to_csv(sys.stdout, sep="\t", index=False)
        return 0
    if not (a.t1 and a.brain and a.lesions):
        ap.error("give t1, brain and lesions, or --manifest")
    df, s = detect(a.t1, a.brain, a.lesions)
    if a.lesions_out:
        df.to_csv(a.lesions_out, sep="\t", index=False)
    print(json.dumps(s, indent=2))
    return 0
