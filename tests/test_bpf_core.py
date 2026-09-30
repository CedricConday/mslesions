import nibabel as nib
import numpy as np

from mslesions.bpf.core import bpf


def test_fraction(tmp_path):
    t1 = np.zeros((20, 20, 20), np.float32)
    t1[2:18, 2:18, 2:18] = 300
    t1[2:18, 2:18, 2:6] = 200
    t1[2:18, 2:18, 14:18] = 50
    les = np.zeros_like(t1, np.uint8)
    les[8:10, 8:10, 8:10] = 1
    t1[les > 0] = 60     # T1-dark lesion: must count as tissue, not CSF
    b = (t1 > 0).astype(np.uint8)
    tis = np.where(b > 0, 3, 0).astype(np.uint8)
    tis[2:18, 2:18, 2:6] = 2
    tis[2:18, 2:18, 14:18] = 1
    for n, x in (("t", t1), ("b", b), ("l", les), ("tis", tis)):
        nib.save(nib.Nifti1Image(x, np.eye(4)), tmp_path / f"{n}.nii.gz")
    r = bpf(tmp_path / "t.nii.gz", tmp_path / "b.nii.gz", tmp_path / "l.nii.gz", tissue_path=tmp_path / "tis.nii.gz")
    assert abs(r["bpf"] - 0.75) < 1e-6 and r["lesion_ml"] > 0
