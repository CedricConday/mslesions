"""mcdis T1w.nii.gz brainmask.nii.gz lesions.nii.gz [--space mni|native] [--it-mask m.nii.gz] | mcdis --manifest m.tsv"""

import argparse
import json
import sys

import pandas as pd

from .core import label


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="mcdis", description=__doc__)
    ap.add_argument("t1", nargs="?")
    ap.add_argument("brain", nargs="?")
    ap.add_argument("lesions", nargs="?")
    ap.add_argument("--space", choices=("mni", "native"), required=True, help="mni uses the MNI ventricle seed and infratentorial box")
    ap.add_argument("--tissue", help="external tissue segmentation on the T1w grid (0 bg, 1 CSF, 2 GM, 3 WM)")
    ap.add_argument("--it-mask", help="infratentorial mask on the T1w grid (native space)")
    ap.add_argument("--manifest", help="TSV: subject, t1, brainmask, mask; writes one row per subject")
    ap.add_argument("--lesions-out", help="per-lesion TSV (single-subject mode)")
    a = ap.parse_args(argv)
    if a.manifest:
        m = pd.read_csv(a.manifest, sep="\t")
        rows = []
        for _, r in m.iterrows():
            try:
                _, s = label(r["t1"], r["brainmask"], r["mask"], space=a.space, it_mask_path=a.it_mask,
                             tissue_path=r.get("tissue") if "tissue" in m else None)
            except ValueError as e:  # TissueError included: a failed tissue split is reported, not guessed
                rows.append({"subject": r["subject"], "error": str(e)})
                print(r["subject"], "FAILED", e, file=sys.stderr, flush=True)
                continue
            rows.append({"subject": r["subject"], **{k: v for k, v in s.items() if k != "note"}})
            print(r["subject"], s["brain_topographies"], file=sys.stderr, flush=True)
        pd.DataFrame(rows).to_csv(sys.stdout, sep="\t", index=False)
        return 0
    if not (a.t1 and a.brain and a.lesions):
        ap.error("give t1, brain and lesions, or --manifest")
    df, s = label(a.t1, a.brain, a.lesions, space=a.space, it_mask_path=a.it_mask, tissue_path=a.tissue)
    if a.lesions_out:
        df.to_csv(a.lesions_out, sep="\t", index=False)
    print(json.dumps(s, indent=2))
    return 0
