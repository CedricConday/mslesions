"""T1 black holes from a T1w, a brain mask and a T2/FLAIR lesion mask.

A lesion is a black hole when the median T1w intensity of its core (voxels at least ``core_mm`` inside
the lesion, or the whole lesion when it is thinner) is at or below the median of cortical-level grey
matter: the "isointense to grey matter or darker" reading used in MAGNIMS guidance and most trials.
Grey matter is the middle class of a three-class k-means on the T1w inside the brain, lesions excluded.
Tissue classes come from the gated split in blackholes.tissue (or an external segmentation); a scan whose
split is implausible raises TissueError instead of returning a threshold. ``persistent`` counts black holes of
one visit that overlap a black hole of a later visit on the same grid.
"""

from __future__ import annotations

import nibabel as nib
import numpy as np
import pandas as pd
from scipy import ndimage as ndi


def _kmeans3(x, iters=30):
    c = np.percentile(x, [15, 50, 85]).astype(float)
    for _ in range(iters):
        lab = np.argmin(np.abs(x[:, None] - c[None]), axis=1)
        new = np.array([x[lab == k].mean() if (lab == k).any() else c[k] for k in range(3)])
        if np.allclose(new, c):
            break
        c = new
    return lab, c


def detect(t1_path, brain_path, lesion_path, core_mm: float = 1.0, min_mm3: float = 14.0, tissue_path=None,
           gate: bool = True) -> tuple[pd.DataFrame, dict]:
    im = nib.load(t1_path)
    t1 = np.asarray(im.dataobj, dtype=np.float32)
    z = tuple(float(v) for v in im.header.get_zooms()[:3])
    vox = float(np.prod(z))
    b = np.asarray(nib.load(brain_path).dataobj) > 0
    les = np.asarray(nib.load(lesion_path).dataobj) > 0.5
    from msdataqc.tissue import classes

    lab = classes(t1, b, les, z, tissue_path, gate=gate)
    gm = float(np.median(t1[lab == 1]))
    wm = float(np.median(t1[lab == 2]))
    depth = ndi.distance_transform_edt(les, sampling=z)
    ll, _n = ndi.label(les, ndi.generate_binary_structure(3, 2))
    rows = []
    for i, sl in enumerate(ndi.find_objects(ll), start=1):
        m = ll[sl] == i
        v = float(m.sum() * vox)
        core = m & (depth[sl] > core_mm)  # the lesion's surface voxels have depth 1 voxel; the core is deeper
        vals = t1[sl][core if core.any() else m]
        med = float(np.median(vals))
        rows.append({"lesion": i, "volume_mm3": v, "counted": v >= min_mm3, "t1_median": med,
                     "t1_to_wm": med / wm, "t1_to_gm": med / gm, "black_hole": bool(v >= min_mm3 and med <= gm)})
    df = pd.DataFrame(rows, columns=["lesion", "volume_mm3", "counted", "t1_median", "t1_to_wm", "t1_to_gm", "black_hole"])
    c = df[df["counted"]]
    bh = c[c["black_hole"]]
    s = {"lesions": len(c), "black_holes": len(bh), "t2_volume_mm3": float(c["volume_mm3"].sum()),
         "black_hole_volume_mm3": float(bh["volume_mm3"].sum()),
         "black_hole_share": float(bh["volume_mm3"].sum() / c["volume_mm3"].sum()) if len(c) else None,
         "gm_median": gm, "wm_median": wm}
    return df, s


def persistent(labels_a, labels_b, bh_a: set, bh_b: set) -> int:
    """Black holes of visit a whose voxels overlap a black hole of visit b (label maps on one grid)."""
    a = np.asarray(nib.load(labels_a).dataobj).astype(np.int64)
    b = np.asarray(nib.load(labels_b).dataobj).astype(np.int64)
    keep = 0
    for i in bh_a:
        over = set(np.unique(b[a == i])) - {0}
        keep += bool(over & bh_b)
    return keep
