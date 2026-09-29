"""Optic nerve as the fifth McDonald 2024 topography.

Rules as recommended by Saidha et al., Lancet Neurol 2025 (the OCT/VEP consensus accompanying the 2024
McDonald criteria): an optic nerve lesion is supported by an OCT inter-eye difference (IED) of at least
6 um in peripapillary RNFL or at least 4 um in macular GCIPL, or by a delayed or asymmetric VEP P100
latency against centre-specific norms (2.5 SD). Conditions the consensus attaches, enforced here:
OCT IEDs are validated from 3 months after acute optic neuritis; scans must pass quality control
(OSCAR-IB); bilateral involvement is not validated by IED; and no better explanation (other eye disease,
high refractive error, uncontrolled diabetes or hypertension) may exist. The last one is the
clinician's; the tool records whether it was confirmed.
"""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

TOPOGRAPHIES = ("optic_nerve", "juxtacortical_cortical", "periventricular", "infratentorial", "spinal_cord")


@dataclass(frozen=True)
class Rules:
    prnfl_ied_um: float = 6.0
    gcipl_ied_um: float = 4.0
    min_days_after_on: int = 90
    vep_sd: float = 2.5


def optic_nerve(oct_row: dict, vep_z: float | None = None, last_on_date=None, exam_date=None,
                other_cause_excluded: bool | None = None, r: Rules | None = None) -> dict:
    """oct_row: prnfl_od, prnfl_os, gcipl_od, gcipl_os (um), qc_pass (bool). Returns status and reasons."""
    r = r or Rules()
    reasons, notes = [], []
    status = None
    if oct_row:
        if oct_row.get("bilateral_disease"):
            notes.append("bilateral optic nerve involvement suspected: inter-eye differences are not validated for it")
        elif oct_row.get("qc_pass") is not True:
            notes.append("OCT quality control not passed or not recorded (OSCAR-IB): not assessable")
        elif last_on_date is not None and exam_date is not None and (pd.Timestamp(exam_date) - pd.Timestamp(last_on_date)).days < r.min_days_after_on:
            notes.append(f"OCT less than {r.min_days_after_on} days after acute optic neuritis: IED not validated (swelling)")
        else:
            compared = False
            for layer, thr in (("prnfl", r.prnfl_ied_um), ("gcipl", r.gcipl_ied_um)):
                a, b = oct_row.get(f"{layer}_od"), oct_row.get(f"{layer}_os")
                if a is None or b is None or pd.isna(a) or pd.isna(b):
                    continue
                compared = True
                ied = abs(a - b)
                if ied >= thr:
                    reasons.append(f"{layer.upper()} inter-eye difference {ied:.1f} um >= {thr:g}")
                    status = True
            if status is None and compared:
                status = False
    if vep_z is not None and not pd.isna(vep_z):
        if vep_z >= r.vep_sd:
            reasons.append(f"VEP latency {vep_z:.1f} SD above centre norms")
            status = True
        elif status is None:
            status = False
    if status and other_cause_excluded is False:
        notes.append("an alternative explanation was not excluded: not counted")
        status = None
    elif status and other_cause_excluded is None:
        notes.append("counted provisionally: confirm no better explanation (eye disease, refraction, diabetes, hypertension)")
    return {"optic_nerve": status, "reasons": reasons, "notes": notes}


def dis(topographies: dict) -> dict:
    """topographies: name -> True/False/None for the five McDonald 2024 locations."""
    pos = [t for t in TOPOGRAPHIES if topographies.get(t) is True]
    unknown = [t for t in TOPOGRAPHIES if topographies.get(t) is None]
    met = len(pos) >= 2
    return {"positive": pos, "not_assessed": unknown, "count": len(pos), "dis_met": met,
            "could_change": (not met) and len(pos) + len(unknown) >= 2}
