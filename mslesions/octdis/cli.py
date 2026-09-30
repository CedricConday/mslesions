"""octdis --prnfl OD OS --gcipl OD OS [--vep-z Z] [--on-date D --exam-date D] [--mcdis summary.json] [--spinal yes|no]"""

import argparse
import json

from .core import dis, optic_nerve


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="octdis", description=__doc__)
    ap.add_argument("--prnfl", nargs=2, type=float, metavar=("OD", "OS"))
    ap.add_argument("--gcipl", nargs=2, type=float, metavar=("OD", "OS"))
    ap.add_argument("--qc-fail", action="store_true")
    ap.add_argument("--bilateral", action="store_true", help="bilateral optic nerve disease suspected: IED not interpretable")
    ap.add_argument("--vep-z", type=float)
    ap.add_argument("--on-date")
    ap.add_argument("--exam-date")
    ap.add_argument("--other-cause-excluded", choices=("yes", "no"))
    ap.add_argument("--mcdis", help="mcdis JSON summary for the brain topographies")
    ap.add_argument("--spinal", choices=("yes", "no"))
    a = ap.parse_args(argv)
    row = {}
    if a.prnfl:
        row.update(prnfl_od=a.prnfl[0], prnfl_os=a.prnfl[1])
    if a.gcipl:
        row.update(gcipl_od=a.gcipl[0], gcipl_os=a.gcipl[1])
    row["qc_pass"] = not a.qc_fail  # the CLI user asserts QC by not passing --qc-fail
    row["bilateral_disease"] = a.bilateral
    oc = None if a.other_cause_excluded is None else a.other_cause_excluded == "yes"
    on = optic_nerve(row, a.vep_z, a.on_date, a.exam_date, oc)
    out = {"optic_nerve": on}
    topo = {"optic_nerve": on["optic_nerve"], "spinal_cord": None if a.spinal is None else a.spinal == "yes"}
    if a.mcdis:
        with open(a.mcdis) as fh:
            s = json.load(fh)
        assessed = s.get("assessed", ["periventricular", "juxtacortical_cortical"] +
                         ([] if "Infratentorial not assessed" in s.get("note", "") else ["infratentorial"]))
        for t in ("juxtacortical_cortical", "periventricular", "infratentorial"):
            topo[t] = (t in s["brain_topographies"]) if t in assessed else None
    out["dis_2024"] = dis(topo)
    print(json.dumps(out, indent=2))
    return 0
