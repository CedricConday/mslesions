import nibabel as nib
import numpy as np

from mslesions.dis.core import Params, label


def _phantom(tmp):
    sh = (100, 120, 100)
    z, y, x = np.indices(sh)
    c = np.array(sh) / 2
    r = np.sqrt(((z - c[0]) / 42) ** 2 + ((y - c[1]) / 52) ** 2 + ((x - c[2]) / 42) ** 2)
    brain = r < 1
    t1 = np.where(brain, 300.0, 0).astype(np.float32)          # white matter
    t1[brain & (r > 0.9)] = 200                                  # cortex
    vent = (np.abs(z - c[0]) < 6) & (np.abs(y - c[1]) < 14) & (np.abs(x - c[2]) < 4)
    t1[vent] = 60                                                # ventricle
    t1[brain & (r > 0.97)] = 60                                  # outer CSF
    les = np.zeros(sh, np.uint8)
    les[50:54, 60:64, 54:58] = 1                                 # abuts the ventricle (x 46..53)
    les[50:54, 60:64, 86:89] = 1                                 # under the cortex
    les[30:34, 60:64, 30:34] = 1                                 # deep white matter
    t1[les > 0] = 150
    tis = np.where(brain, 3, 0).astype(np.uint8)
    tis[brain & (r > 0.9)] = 2
    tis[vent | (brain & (r > 0.97))] = 1
    aff = np.eye(4)
    for n, a in (("t1", t1), ("brain", brain.astype(np.uint8)), ("les", les), ("tis", tis)):
        nib.save(nib.Nifti1Image(a, aff), tmp / f"{n}.nii.gz")
    return tmp / "t1.nii.gz", tmp / "brain.nii.gz", tmp / "les.nii.gz", tmp / "tis.nii.gz"


def test_classes(tmp_path):
    t1, b, les, tis = _phantom(tmp_path)
    df, s = label(t1, b, les, Params(min_lesion_mm3=10.0), space="native", tissue_path=tis)
    assert sorted(df["class"]) == ["juxtacortical_cortical", "other", "periventricular"]
    assert s["dis_2024_from_brain_mri"] and s["n_brain_topographies"] == 2
    assert "Infratentorial not assessed" in s["note"]
