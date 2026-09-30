# mslesions

Lesion measures from MS MRI, each behind a tissue-split gate that refuses a scan it cannot read instead of returning a number (the gate lives in [msdataqc](https://github.com/CedricConday/msdataqc)).

```bash
pip install mslesions
mslesions <tool> --help
```

| Tool | Command | What it does |
|---|---|---|
| [`dis`](#dis) | `mslesions dis` | McDonald 2024 topography of every lesion and the DIS count (formerly dis) |
| [`octdis`](#octdis) | `mslesions octdis` | optic nerve by OCT/VEP and the five-location DIS count |
| [`blackholes`](#blackholes) | `mslesions blackholes` | T1 black holes from a T2 lesion mask |
| [`bpf`](#bpf) | `mslesions bpf` | brain parenchymal fraction with an age self-test (experimental) |

Research software, not a medical device. MIT licence. Each section below keeps the tool's own results and the
corrections that three independent reviews made to them.

## dis

*formerly mcdis*

**Dissemination in space from a lesion mask.** dis labels every multiple sclerosis lesion as
periventricular, juxtacortical/cortical, infratentorial or other, and counts the McDonald 2024
topographies a brain MRI shows. Input: a T1w, a brain mask and a lesion mask (from an expert or a
segmenter such as LST-AI). No atlas registration for the supratentorial classes; infratentorial uses
an MNI-space rule or a mask you supply.

```bash
pip install mslesions
mslesions dis T1w.nii.gz brainmask.nii.gz lesions.nii.gz --lesions-out lesions.tsv     # MNI space
mslesions dis T1w.nii.gz brain.nii.gz lesions.nii.gz --space native --it-mask infratentorial.nii.gz
mslesions dis --manifest cohort.tsv > dis.tsv                                          # subject, t1, brainmask, mask
```

How: three-class k-means on the T1w (lesion voxels counted as white matter); ventricles are CSF more
than 12 mm inside the brain; cortex is grey matter within 6 mm of the surface; a lesion within 2 mm of
a ventricle is periventricular, else within 1.5 mm of cortex juxtacortical/cortical, else
infratentorial if its centre lies below the tentorium. Lesions under 10 mm³ are listed but not counted.
All thresholds are in `mslesions/dis/core.py`.

### MSLesSeg baselines (corrected 2026-09-30)

An earlier version reported 36 of 53 patients with dissemination in space. An independent review found three
faults behind it: the tissue split put extracranial fat into "white matter", the "ventricles" were all deep CSF
(fissures and cisterns included), and the infratentorial box reached into the fusiform and lingual gyri while
brainstem lesions were taken by the periventricular rule first. All three are fixed: the tissue split is gated
(refuses implausible classes, or takes an external segmentation), the ventricles are the CSF connected to a seed
on the lateral ventricle bodies, infratentorial is tested first with a conservative box that stays below the
temporal and occipital cortex, and lesions under 14 mm³ (a 3 mm sphere) are not counted.

| | Patients |
|---|---:|
| measurable | 28 of 53 |
| refused: grey and white matter do not separate on the T1w | 23 |
| refused: ventricle estimate implausible | 2 |
| periventricular, juxtacortical/cortical, infratentorial | 11 |
| periventricular, infratentorial | 7 |
| periventricular, juxtacortical/cortical | 5 |
| periventricular only | 5 |

23 of the 28 fulfil dissemination in space from brain MRI alone; the other 5 would need the optic nerve
([octdis](https://github.com/CedricConday/mslesions#octdis)) or the spinal cord. Ventricle volume median 63 ml in MNI space.
The infratentorial box is conservative: it misses the upper vermis and the superior cerebellar surface, so an
atlas-derived mask (`--it-mask`) is the right input for a real analysis. No expert location labels exist for
MSLesSeg to measure the rest.

## octdis


McDonald 2024 made the **optic nerve the fifth topography** for dissemination in space, and allowed OCT and
VEP to show it. octdis applies the consensus rules (Saidha et al., Lancet Neurol 2025): an inter-eye
difference of at least 6 µm in peripapillary RNFL or at least 4 µm in macular GCIPL, or a VEP latency
2.5 SD beyond centre norms. It enforces the consensus conditions: IEDs only from 3 months after acute optic
neuritis, quality-controlled scans, and no better explanation. Combined with brain-MRI topographies from
[dis](https://github.com/CedricConday/mslesions#dis) and a spinal-cord answer, it returns the five-location
count, the locations not yet assessed, and whether assessing them could still change the answer.

```bash
pip install mslesions
mslesions octdis --prnfl 96 88 --gcipl 71 69 --exam-date 2025-06-01 --on-date 2025-01-10 --mcdis P1_dis.json --spinal no
```

Of the 17 MSLesSeg patients whom dis finds with only one brain topography, every one is `could_change`.
An OCT is the cheapest test that could complete their dissemination in space. That is the practical point.

## blackholes


**T1 black holes** from a T1w, a brain mask and a T2/FLAIR lesion mask. A lesion is a black hole when
the median T1w intensity of its core is at or below the grey-matter median of the same scan. That is the
"isointense to grey matter or darker" reading used in trials, and it needs no absolute threshold. Output:
per-lesion T1 ratios, black-hole count, volume and share, plus persistence across visits.

```bash
pip install mslesions
mslesions blackholes T1w.nii.gz brainmask.nii.gz lesions.nii.gz --lesions-out lesions.tsv
mslesions blackholes --manifest cohort.tsv > burden.tsv        # subject, t1, brainmask, mask
```

### MSLesSeg baselines (corrected 2026-09-30)

An earlier version reported black-hole shares near 0.97 for some patients. An independent review found why: the
tissue split behind the grey-matter threshold put extracranial fat and venous sinus into its brightest class, so
the "grey matter" median was really the median of all brain tissue. The split now fits on the brain eroded by
5 mm and refuses implausible classes, and the lesion core really is the interior (deeper than one voxel).

With that gate, 30 of the 53 baselines can be measured; the other 23 have T1w scans in which grey and white
matter do not separate (see [intensitydrift](https://github.com/CedricConday/msdataqc#intensitydrift)). On the 30:

| | Value |
|---|---:|
| black holes per patient (median) | 4.5 |
| share of T2 lesion volume, median per patient | 25 % |
| share of T2 lesion volume, pooled | 52 % |
| Spearman with baseline EDSS (0 to 6.0): T2 volume | 0.20 (p 0.30) |
| Spearman with baseline EDSS: black-hole volume | 0.19 (p 0.30) |

Neither burden tracks EDSS in these 30 patients, which is too few to say more. The pooled share is high because a
few large confluent lesions are dark throughout. Give an external tissue segmentation when the built-in split is
refused.

## bpf


Brain parenchymal fraction, tissue over intracranial volume, from a T1w and a brain mask, with T1-dark
lesions counted as tissue. It needs no atlas and no external segmenter.

```bash
pip install mslesions
mslesions bpf T1w.nii.gz brainmask.nii.gz lesions.nii.gz
mslesions bpf --manifest cohort.tsv > bpf.tsv          # subject, t1, brainmask[, mask][, age]: prints the age check
```

### The check it failed

BPF falls with age in every healthy and MS cohort, so the manifest mode prints the correlation with age as a
self-test. On the 53 MSLesSeg baselines it is **-0.07**: no relation. The per-class volumes show why. The
three-class clustering of the T1w gives white matter anywhere from 161 to 1004 ml across patients, because
MSLesSeg's T1w images are MNI-normalised and their intensity distributions differ a lot between patients.
Clipping the bright tail (vessels, fat) helped, but not enough. Since 2026-09-30 bpf uses the gated tissue split
shared with dis, intensitydrift and blackholes: 23 of the 53 scans are refused (their T1w does not separate
grey and white matter), and on the 30 that pass the correlation with age is 0.00. The self-test still fails.

So on data like MSLesSeg this tool does not measure atrophy, and says so. Use it on native-space T1w from one
scanner, or replace the clustering with a real tissue segmenter (FSL FAST, SynthSeg). It stays here because the
failed self-test is the useful part: a cross-sectional atrophy number without an age check is not to be trusted.
