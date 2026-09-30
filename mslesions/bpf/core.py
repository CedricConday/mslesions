"""Brain parenchymal fraction (BPF) = (grey + white matter + lesions) / (grey + white + lesions + CSF),
inside the brain mask, from a three-class k-means on the T1w clipped to its 1st-95th percentiles (lesion voxels counted as tissue,
since T1-dark lesions would otherwise read as CSF; the clip keeps bright vessels from taking a class). The brain mask defines the intracranial space, so a mask that
excludes outer CSF gives a BPF that is too high; the ``mask_note`` field says which kind was used.
"""

from __future__ import annotations

import nibabel as nib
import numpy as np


def _kmeans3(x, iters=40):
    c = np.percentile(x, [15, 50, 85]).astype(float)
    for _ in range(iters):
        lab = np.argmin(np.abs(x[:, None] - c[None]), axis=1)
        new = np.array([x[lab == k].mean() if (lab == k).any() else c[k] for k in range(3)])
        if np.allclose(new, c):
            break
        c = new
    return lab


def bpf(t1_path, brain_path, lesion_path=None, tissue_path=None, gate: bool = True) -> dict:
    im = nib.load(t1_path)
    t1 = np.asarray(im.dataobj, dtype=np.float32)
    b = np.asarray(nib.load(brain_path).dataobj) > 0
    les = np.asarray(nib.load(lesion_path).dataobj) > 0.5 if lesion_path else np.zeros_like(b)
    vox = float(np.prod(im.header.get_zooms()[:3])) / 1000.0
    from msdataqc.tissue import classes

    lab = classes(t1, b, les, tuple(float(z) for z in im.header.get_zooms()[:3]), tissue_path, gate=gate)
    csf = float((lab == 0).sum() * vox)
    gm = float((lab == 1).sum() * vox)
    wm = float((lab == 2).sum() * vox)
    lv = float((les & b).sum() * vox)
    icv = csf + gm + wm + lv
    return {"bpf": (gm + wm + lv) / icv, "gm_ml": gm, "wm_ml": wm, "lesion_ml": lv, "csf_ml": csf, "icv_ml": icv,
            "gm_fraction": gm / icv, "wm_fraction": (wm + lv) / icv}
