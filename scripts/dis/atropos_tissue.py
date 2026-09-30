"""External tissue segmentation for mcdis / intensitydrift / blackholes: ANTs Atropos (3 classes, k-means
initialisation, MRF smoothing) inside the brain mask with lesions excluded, after N4 bias correction.
Writes 0 background, 1 CSF, 2 GM, 3 WM on the T1w grid. Usage: atropos_tissue.py T1 brainmask lesions out.nii.gz
Needs antspyx (pip install antspyx)."""
import sys

import ants
import numpy as np

t1, bm, les, out = sys.argv[1:5]
img = ants.image_read(t1)
brain = ants.image_read(bm)
lesion = ants.image_read(les)
mask_np = (brain.numpy() > 0) & (lesion.numpy() <= 0.5)
mask = brain.new_image_like(mask_np.astype("float32"))
n4 = ants.n4_bias_field_correction(img, mask=mask)
seg = ants.atropos(a=n4, x=mask, i="kmeans[3]", m="[0.2,1x1x1]", c="[5,0]")
lab = seg["segmentation"].numpy().astype(np.uint8)  # 1..3 ordered by mean intensity: CSF, GM, WM on T1w
ants.image_write(brain.new_image_like(lab.astype("float32")), out)
v = np.prod(img.spacing) / 1000
print({k: round(float((lab == k).sum() * v)) for k in (1, 2, 3)})
