# eclipse-hdr — standalone program

Solar-eclipse HDR processing **without Siril**: all your RAW files in one
folder, a pipeline of independent steps, a linear 32-bit TIFF at every step, and
a 16-bit export at any point for Photoshop, DxO or Affinity.

> **Status: phase 1.** Inventory, decoding, stacking and export work and are
> validated on real RAW files. Alignment, HDR merge, corona, compositing and the
> graphical interface still have to be ported from the Siril scripts — see
> [What is still missing](#what-is-still-missing).

---

## Installation

```bash
cd Eclipse_HDR
python3 -m venv .venv
.venv/bin/pip install -e .
```

That installs everything: numpy, scipy, opencv, rawpy, tifffile, astropy,
exifread. Nothing else is needed — no Siril, no Photoshop, not even exiftool
(which only helps read EXIF from recent containers such as CR3:
`brew install exiftool`).

## Quick start

```bash
.venv/bin/eclipse-hdr scan   ~/totality
.venv/bin/eclipse-hdr decode ~/totality
.venv/bin/eclipse-hdr stack  ~/totality
.venv/bin/eclipse-hdr status ~/totality
```

`~/totality` is the **project folder**. It holds the RAW files, either directly
or in a `raw/` subfolder. Nothing is written anywhere else, and the RAW files
are never modified, moved or copied.

## The project folder

```
totality/
  raw/                  ← your RAW files (or directly in totality/)
  projet.json           ← manifest: groups, parameters, state of each step
  01_decode/1-160/…     ← linear + saturation mask, sorted by shutter speed
  03_stack/master_1-160.tif
  export/
```

The manifest is the **logbook** of the run: which speeds were found, what
settings each step used, when, and how long it took. It is what makes resuming
possible.

**Incremental resume.** Each step records a fingerprint
`sha256(name + parameters + inputs)`. Re-running the pipeline only recomputes
what has changed: tweaking a corona setting re-stacks nothing. Files already
written are skipped.

**Interruption.** Ctrl-C stops after the current file (a second Ctrl-C forces).
Every write is atomic: an interruption never leaves a truncated file that
*looks* valid.

---

## The steps

| Command | Purpose | Output |
|---|---|---|
| `scan` | read EXIF, group by shutter speed | manifest |
| `decode` | RAW → linear + saturation mask | `01_decode/<speed>/` |
| `stack` | one master per speed | `03_stack/master_<speed>.tif` |
| `export` | 16-bit TIFF for Photoshop/DxO | anywhere you like |
| `info` | metadata of an image | terminal |
| `diff` | compare two images | terminal |
| `status` | project state | terminal |
| `run` | chain the available steps | — |

### `scan` — inventory

```bash
eclipse-hdr scan ~/totality --tol 3
```

Reads each RAW file's exposure time **without converting anything**: a few
seconds, safe to re-run to inspect the plan before committing. `--tol` is the
grouping tolerance in percent (3 % by default: merges APEX rounding such as
1/1000 vs 1/1024, without mixing two ⅓-stop steps, 26 % apart).

Speeds become folder names: 1/160 s → `1-160`, 2 s → `2s`, 1.3 s → `1s3`.

### `decode` — RAW to linear

```bash
eclipse-hdr decode ~/totality --speeds 1-160,2s --jobs 8
```

| Option | Default | Effect |
|---|---|---|
| `--speeds` | all | comma-separated speeds to process |
| `--wb` | `daylight` | white balance, the same for every image |
| `--jobs` | cores (max 8) | parallel processes |
| `--demosaic` | `AHD` | `AHD`, `VNG`, `PPG`, `DHT`, `AAHD`, `LINEAR` |
| `--sat-frac` | `0.98` | saturation threshold, as a fraction of the usable range |
| `--sat-dilate` | `3` | mask dilation in px (demosaic reach) |
| `--half` | off | quarter resolution, for a quick trial |
| `--f32` | off | float32 instead of 16-bit (see below) |
| `--force` | off | redo even if the file exists |

Budget roughly **1.5 s per image per core**. Four hundred images on eight cores:
about a minute and a half.

#### `--wb` — fixing the white balance

A **single** balance is applied to every image of the session, and it is saved
in the manifest: later decode runs reuse it without you having to repeat it.

| Value | Balance used |
|---|---|
| *(omitted)* or `daylight` | the camera's daylight balance — **default** |
| `camera` | as shot, read from the first file |
| `neutral` | none: equal multipliers |
| `1.9722,0.9412,1.1376` | explicit values, three (RGB) or four (RGBG) |
| `~/photos/_DSC0577.NEF` | the as-shot balance of **that** file |

To see what a RAW file holds before choosing:

```bash
eclipse-hdr info ~/totality/raw/_DSC0577.NEF
```

```
pose : 1/640 s  (0.0015625 s)
balance lumière du jour  : 1.9722,0.9412,1.1376,0.9412   → --wb 1.9722,0.9412,1.1376
balance telle que prise  : 1.7930,1.0000,1.5352,1.0000   → --wb 1.7930,1.0000,1.5352
```

> **Only the ratios matter.** libraw divides the four multipliers by the
> smallest before applying them: `2,2,2,2` yields exactly the same image as
> `1,1,1,1` (verified). This setting changes the **tint**, never the brightness.

The default is the daylight balance because it is a **constant of the camera
body** — see [The radiometric invariant](#the-radiometric-invariant) for the
measurement that makes it necessary.

**Changing your mind later.** If images are already decoded, a different `--wb`
is **refused**:

```
✗ la balance des blancs change (…) alors que 2 vitesse(s) sont déjà décodées.
  Mélanger deux balances dans un même projet casse l'échelle radiométrique
  dont dépend la fusion HDR.
  • pour tout redécoder avec la nouvelle balance : ajoutez --force
  • pour garder l'ancienne : relancez sans --wb
```

With `--force`, **every** already-decoded speed is redone, including those
`--speeds` did not ask for: a project can therefore never end up holding two
mixed balances. A balance merely proportional to the previous one is not a
conflict, since it produces the same image.

### `stack` — one master per speed

```bash
eclipse-hdr stack ~/totality
```

Without `--method`, each group gets the method suited to its frame count, and
the choice is logged:

| Frames | Method | Siril equivalent |
|---|---|---|
| 2 – 3 | `median` | `stack med -nonorm` |
| 4 – 6 | `percentile` 0.2/0.1 | `stack rej p 0.2 0.1 -nonorm` |
| 7 and up | `winsorized` 3/3 | `stack rej w 3 3 -nonorm` |

| Option | Default | Effect |
|---|---|---|
| `--method` | auto | `median`, `mean`, `sum`, `min`, `max`, `percentile`, `sigma`, `winsorized` |
| `--low` / `--high` | per method | rejection thresholds |
| `--missing` | `nan` | what counts as missing data (see below) |
| `--band` | `256` | memory budget per band, in MB |

**No normalisation is applied, and there is no way to ask for any**: the
equivalent of `-nonorm` is the only possible behaviour. The trap that breaks the
HDR merge cannot spring.

**Bounded memory.** Images are processed in horizontal bands: peak memory
depends on `--band`, not on the number of frames. Forty 45 Mpx images fit in
300 MB instead of 21 GB.

### `export` — 16-bit for Photoshop and DxO

```bash
eclipse-hdr export 03_stack/master_1-160.tif -o for_photoshop.tif
eclipse-hdr export 05_hdr/hdr.tif --transform asinh --asinh-k 300
```

| `--transform` | When |
|---|---|
| `autostretch` *(default)* | everyday use — colours preserved (the stretch acts on luminance) |
| `asinh` | gentle stretch tunable with `--asinh-k`, keeps the core from blowing out |
| `linear` | scaling only, if you stretch elsewhere |

> ⚠ The default is deliberately **not** `linear`. A linear HDR spans five
> decades; taken to 16-bit without a stretch, the outer corona lands on a few
> dozen levels and the file opens as a black rectangle.

---

## Design decisions, and why

The three claims below were **measured** on real RAW files, not assumed. The
figures come from a batch of twenty NEF files.

### The radiometric invariant

The whole HDR merge rests on *value ∝ exposure × radiance*. Decoding guarantees
it through three settings: `gamma=(1,1)`, `no_auto_bright=True`, and white
balance multipliers **frozen for the whole session**.

The balance used is the camera's daylight one, a constant of the body.
Measurement: across twenty files it takes **one** value, while the *as-shot*
balance takes **twenty different ones**. Using the latter would give every frame
its own scale factor.

Direct check: halving the sensor signal divides the output by **0.4997** (0.03 %
error), with a spread of ±0.0004 over 24.6 million pixels. With auto-brightness
enabled the same test gives **×1.0006** — the exposure difference is erased
entirely.

### The saturation mask, read at the sensor

After white balance and the colour matrix, a saturated photosite no longer lands
on a single value: each channel is multiplied differently, then mixed.
Measurement: pixels saturated at the sensor come out spread over **0.254** in
the output. No threshold on the final image can fence them off cleanly.

The mask is read **before demosaicing**, by comparison with the sensor's
saturation level, then dilated by three pixels — a saturated photosite
contaminates its interpolated neighbours (measured: 277 → 1855 pixels).

The mask travels with the image all the way to the merge, and stacking writes
the **fraction of saturated frames** at every pixel.

### 16-bit storage for the decode step

libraw's output *is* 16-bit integer. Storing it as float32 pads it with zeros:
same information, twice the disk. Over a session of four hundred 45 Mpx images
that is **46 GB instead of 92**. The round trip is exact to the bit (verified).

Later steps produce float32: from stacking onwards, averaging several frames
creates genuine intermediate values. `--f32` makes it uniform if you prefer.

---

## Validation against Siril

The same four images (cropped to 1200×800, i.e. 2,880,000 values) stacked by
Siril 1.4.4 and by `eclipse-hdr`, from input data **identical to the bit**:

| Method | Divergent pixels | Largest difference |
|---|---|---|
| `median` | **0** | 0 — exact to the bit |
| `mean` | **0** | 3·10⁻⁸ (float32 rounding) |
| `percentile` 0.2/0.1 | **1 in 2,880,000** | one value 2·10⁻⁸ from the rejection threshold |
| `winsorized` 3/3 | divergent | see below |

To reproduce: `eclipse-hdr diff my_master.tif siril_master.fit`.

### Two known divergences

**Zero pixels.** Siril excludes exactly-zero pixels from its **mean** stacks —
but not from its **median**. On `[0, 0.2179, 0.2313, 0.2459]` its mean returns
0.2317 (the mean of the three non-zero values) where the exact mean is 0.1738.
That affected 4.3 % of the pixels in the test.

The convention is defensible for registered frames, where the borders left empty
by the shift genuinely hold no data. It is not, for raw frames, where zero is a
legitimate measurement (a photosite below the black level). Here, missing data
is therefore marked **NaN**, never zero, and `--missing zero` reproduces Siril's
convention — that is the setting under which the mean matches exactly.

**Winsorized clipping.** On four frames, Siril's `stack rej w 3 3` rejects
*nothing*: its output is bit-for-bit identical to its plain mean. The one here
does reject, because the winsorized sigma does not let itself be inflated by the
outlier — which is the entire point of it.

The limit can be written down. A value at distance *d* among *n* frames raises
the ordinary standard deviation to √(σ₀² + *d*²/*n*), and is rejected at *k* σ
only if *d*²(1 − *k*²/*n*) > *k*²σ₀². With the usual *k* = 3 the factor vanishes
as soon as **n ≤ 9**: below ten frames, ordinary sigma clipping cannot reject
*any* isolated value, however aberrant. That is what winsorization exists for,
and why the rule says "percentile below seven frames".

This divergence still needs **checking on eclipse frames**: four images of
different scenes are not a representative test case.

---

## Interoperability with Siril

Neither tool forces out the other: the two chains cross.

```bash
eclipse-hdr info  master.fit          # reads FITS too
eclipse-hdr diff  a.tif b.fit         # compares TIFF against FITS
```

In Python, `imageio.from_fits` and `imageio.to_fits` import and export 32-bit
FITS carrying `EXPTIME` and `MOONX`/`MOONY`/`MOONR` — enough to inject Siril
masters at any step, or to go back into Siril for the GHS stretch.

## What is still missing

| Step | State |
|---|---|
| Inventory, decoding, stacking, export | **done, validated on real RAW** |
| Moon (☾) and multi-point (☼) alignment | to port from `SirilJ_Align.py` |
| HDR merge | to port from `FusionHDR_v1.2.py` |
| Corona (RHEF / FNRGF / MGN / ACHF) | to port from `Corona_v1.1.py` |
| Layer compositing | to port from `EclipseComposite.py` |
| PyQt6 graphical interface | to write |

The algorithms already exist and are independent of Siril: porting means
replacing the input/output layer, not rewriting the science.

Pre-processing (bias, darks, flats) is out of scope at this stage, the current
Siril chain not using it either.

## Licence

[GPL-3.0-or-later](../LICENSE).
