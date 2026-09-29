import nibabel as nib
import numpy as np

from mslesions.blackholes.core import detect


def test_dark_lesion_is_black_hole(tmp_path):
    sh = (40, 40, 40)
    t1 = np.zeros(sh, np.float32)
    t1[4:36, 4:36, 4:36] = 300
    t1[4:36, 4:36, 4:10] = 200     # grey matter
    t1[4:36, 4:36, 30:36] = 60     # csf
    les = np.zeros(sh, np.uint8)
    les[14:19, 14:19, 14:19] = 1   # dark: 120 < 200
    les[22:27, 22:27, 14:19] = 1   # bright-ish: 260 > 200
    t1[14:19, 14:19, 14:19] = 120
    t1[22:27, 22:27, 14:19] = 260
    b = (t1 > 0).astype(np.uint8)
    tis = np.where(b > 0, 3, 0).astype(np.uint8)
    tis[4:36, 4:36, 4:10] = 2
    tis[4:36, 4:36, 30:36] = 1
    for n, x in (("t1", t1), ("b", b), ("l", les), ("tis", tis)):
        nib.save(nib.Nifti1Image(x, np.eye(4)), tmp_path / f"{n}.nii.gz")
    _df, s = detect(tmp_path / "t1.nii.gz", tmp_path / "b.nii.gz", tmp_path / "l.nii.gz", tissue_path=tmp_path / "tis.nii.gz")
    assert s["lesions"] == 2 and s["black_holes"] == 1 and abs(s["black_hole_share"] - 0.5) < 1e-9
