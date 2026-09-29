"""Label MS lesions by McDonald topography from a T1w, a brain mask and a lesion mask.

Tissue: three-class k-means on T1w inside the brain (lesion voxels are counted as white matter,
since they are T1-dark and would otherwise read as CSF or grey matter).
Ventricles: CSF deeper than ``ventricle_depth_mm`` from the brain surface, components over
``ventricle_min_ml``. Cortex: grey matter within ``cortex_depth_mm`` of the brain surface.
Lesion classes, first match wins:
* periventricular: within ``pv_mm`` of a ventricle (MAGNIMS: lesion abutting the ventricle);
* juxtacortical/cortical: within ``jc_mm`` of cortex;
* infratentorial: centre below the tentorium, from an MNI-space rule or a user mask;
* other (deep white matter).
McDonald 2024 counts five topographies (optic nerve, juxtacortical/cortical, periventricular,
infratentorial, spinal cord); a brain MRI can show three of them, and the count says so.
"""

from __future__ import annotations

from dataclasses import dataclass

import nibabel as nib
import numpy as np
import pandas as pd
from scipy import ndimage as ndi


@dataclass(frozen=True)
class Params:
    ventricle_depth_mm: float = 15.0
    ventricle_min_ml: float = 0.5
    cortex_depth_mm: float = 6.0
    pv_mm: float = 2.0
    jc_mm: float = 1.5
    min_lesion_mm3: float = 14.0  # MAGNIMS: lesions of at least 3 mm diameter (a 3 mm sphere is 14 mm3)
    connectivity: int = 2
    tissue_gate: bool = True  # refuse implausible tissue splits (turn off only for phantoms)


def kmeans3(x: np.ndarray, iters: int = 30) -> np.ndarray:
    c = np.percentile(x, [15, 50, 85]).astype(np.float64)
    for _ in range(iters):
        lab = np.argmin(np.abs(x[:, None] - c[None, :]), axis=1)
        new = np.array([x[lab == k].mean() if (lab == k).any() else c[k] for k in range(3)])
        if np.allclose(new, c):
            break
        c = new
    return lab


def _xyz(shape: tuple, affine: np.ndarray):
    ijk = np.indices(shape).reshape(3, -1)
    return (affine[:3, :3] @ ijk + affine[:3, 3:4]).reshape((3, *shape))


def infratentorial_mni(shape: tuple, affine: np.ndarray) -> np.ndarray:
    """A conservative cerebellum-and-brainstem box in MNI152 mm. It stays below the inferior temporal and occipital
    cortex (fusiform, lingual and inferior temporal gyri reach z of about -20 lateral to |x| 25), so it misses the
    upper vermis and the superior cerebellar surface rather than calling cortex infratentorial. For a real
    analysis pass an atlas-derived mask with ``it_mask_path``."""
    x, y, z = _xyz(shape, affine)
    brainstem = (np.abs(x) < 12) & (y > -42) & (y < -14) & (z < -12)
    vermis = (np.abs(x) < 10) & (y < -42) & (z < -12)
    hemispheres = (np.abs(x) < 50) & (y < -42) & (z < -26)
    return brainstem | vermis | hemispheres


def lateral_ventricle_seed(shape: tuple, affine: np.ndarray) -> np.ndarray:
    x, y, z = _xyz(shape, affine)
    # bodies of the lateral ventricles, off the midline (the interhemispheric fissure is CSF too)
    return (np.abs(x) > 6) & (np.abs(x) < 24) & (y > -30) & (y < 15) & (z > 12) & (z < 30)


def label(t1_path, brain_path, lesion_path, p: Params | None = None, space: str = "mni", it_mask_path=None,
          tissue_path=None) -> tuple[pd.DataFrame, dict]:
    p = p or Params()
    t1i = nib.load(t1_path)
    zooms = tuple(float(z) for z in t1i.header.get_zooms()[:3])
    vox = float(np.prod(zooms))
    t1 = np.asarray(t1i.dataobj, dtype=np.float32)
    brain = np.asarray(nib.load(brain_path).dataobj) > 0
    les = np.asarray(nib.load(lesion_path).dataobj) > 0.5
    for q in (brain, les):
        if q.shape != t1.shape:
            raise ValueError("T1w, brain mask and lesion mask must share one grid")
    from msdataqc.tissue import classes

    lab = classes(t1, brain, les, zooms, tissue_path, gate=p.tissue_gate)
    lab[les] = 2
    csf, gm = lab == 0, lab == 1
    depth = ndi.distance_transform_edt(brain, sampling=zooms)
    vent_c = csf & (depth > p.ventricle_depth_mm)
    cc, n = ndi.label(vent_c)
    if space == "mni":
        # the lateral ventricles: CSF components that reach the ventricle seed region; the fissures and cisterns do not
        # the ventricles are the CSF under the seed, cut from the rest by opening (thin sulcal bridges break)
        vent_c = ndi.binary_opening(vent_c, iterations=1)
        cc, n = ndi.label(vent_c)
        seed_ids = set(np.unique(cc[lateral_ventricle_seed(t1.shape, t1i.affine)])) - {0}
        sizes = np.bincount(cc.ravel()) * vox / 1000.0
        vent = np.isin(cc, [i for i in seed_ids if sizes[i] >= p.ventricle_min_ml])
    else:
        sizes = np.bincount(cc.ravel()) * vox / 1000.0
        vent = cc == (int(np.argmax(sizes[1:])) + 1) if n else np.zeros_like(csf)
    if vent.sum() * vox / 1000.0 > 150:
        raise ValueError(f"ventricle estimate {vent.sum() * vox / 1000.0:.0f} ml: implausible (sulcal CSF connected?)")
    cortex = gm & (depth <= p.cortex_depth_mm)
    d_vent = ndi.distance_transform_edt(~vent, sampling=zooms) if vent.any() else np.full(t1.shape, np.inf)
    d_ctx = ndi.distance_transform_edt(~cortex, sampling=zooms)
    if it_mask_path is not None:
        it = np.asarray(nib.load(it_mask_path).dataobj) > 0
    elif space == "mni":
        it = infratentorial_mni(t1.shape, t1i.affine)
    else:
        it = None
    ll, nl = ndi.label(les, ndi.generate_binary_structure(3, p.connectivity))
    rows = []
    for i, sl in enumerate(ndi.find_objects(ll), start=1):
        m = ll[sl] == i
        vol = float(m.sum() * vox)
        dv, dc = float(d_vent[sl][m].min()), float(d_ctx[sl][m].min())
        cen = [float(c) for c in ndi.center_of_mass(m)]
        cen = tuple(round(c + s.start) for c, s in zip(cen, sl))
        if it is not None and it[cen]:
            cls = "infratentorial"  # first: a brainstem lesion next to the fourth ventricle is infratentorial
        elif dv <= p.pv_mm:
            cls = "periventricular"
        elif dc <= p.jc_mm:
            cls = "juxtacortical_cortical"
        else:
            cls = "other"
        rows.append({"lesion": i, "volume_mm3": vol, "class": cls, "counted": vol >= p.min_lesion_mm3,
                     "dist_ventricle_mm": dv, "dist_cortex_mm": dc, "centre_vox": cen})
    df = pd.DataFrame(rows, columns=["lesion", "volume_mm3", "class", "counted", "dist_ventricle_mm", "dist_cortex_mm", "centre_vox"])
    counted = df[df["counted"]]
    tops = sorted(set(counted["class"]) - {"other"})
    summary = {
        "lesions": int(nl), "counted": len(counted),
        "by_class": counted["class"].value_counts().to_dict(),
        "brain_topographies": tops, "n_brain_topographies": len(tops),
        "assessed": ["periventricular", "juxtacortical_cortical"] + (["infratentorial"] if it is not None else []),
        "dis_2024_from_brain_mri": len(tops) >= 2,
        "note": "McDonald 2024 DIS: >= 1 typical lesion in >= 2 of 5 topographies (optic nerve, juxtacortical/cortical, "
                "periventricular, infratentorial, spinal cord). A brain MRI shows three; optic nerve and spinal cord "
                "need their own imaging and would only add to this count." + ("" if it is not None else
                " Infratentorial not assessed: native space without an infratentorial mask."),
        "ventricle_ml": float(vent.sum() * vox / 1000.0), "cortex_ml": float(cortex.sum() * vox / 1000.0),
    }
    return df, summary
