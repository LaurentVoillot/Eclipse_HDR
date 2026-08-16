# SirilJ Align — Sun & Eclipse alignment

**One tool, two targets.** A Python script for Siril (`scripts/SirilJ_Align.py`, v1.0.0)
that unifies and replaces two earlier tools: multi-point alignment of the solar surface,
and registration of eclipse exposures on the Moon.

| Mode | Use case | Reference feature |
|---|---|---|
| **☼ Sun — surface** | constant-exposure series (H-alpha, white light, granulation) | surface structure, with local seeing correction |
| **☾ Eclipse — Moon** | exposures of very different brightness (diamond ring to full corona) | the lunar disc, the only feature common to all levels |

In both modes the tool loads **the sequence currently open in Siril** by default.

> The user interface is in French. This guide gives the English meaning of every
> control, with the on-screen label in *italics*.

---

## Why two modes

The two problems are opposites, hence two algorithms.

**Sun.** All frames look alike, but seeing distorts them *locally*: a single global
transform (shift, rotation) cannot correct a distortion that varies across the frame.
The shift has to be measured at many points.

**Eclipse.** Frames do not look alike at all: between a diamond-ring exposure and a
full-corona exposure the content differs too much for Siril's usual correlation to work.
The **lunar limb** is the only feature present and identical in every frame.

---

## Requirements

- Siril 1.3 or newer with Python script support (`sirilpy`).
- Dependencies installed automatically: `numpy`, `astropy`, `scipy`, `opencv-python`, `PyQt6`.
- Frames must be **individual FITS files**. A single-file sequence (FITSEQ/SER) is detected
  and reported: export it to individual FITS first.

---

## ☼ Sun mode — surface (multi-point)

### How it works

1. **Global shift** measured by phase correlation (sub-pixel accuracy).
2. A **grid of anchor points** is placed on genuinely structured areas of the reference
   frame — granulation, prominences, limb. Empty sky is excluded automatically (a point is
   kept only if its content rises above the noise).
3. **Local shift** measured at every anchor point, then a dense, smoothed deformation
   field, with outlier points rejected and filled in from their neighbours.
4. **Refinement passes**: each pass measures the residual and **composes the maps** — so
   the image is resampled **only once**, with no cumulative blurring.

### Steps

1. Open the sequence in Siril, then launch the script (it detects the sequence; otherwise
   *⟳ Rafraîchir* — Refresh).
2. Choose the **reference frame** (*Image de référence*): *Meilleure netteté* (sharpest,
   default), *Première* (first), *Dernière* (last) or *Manuel*.
3. *Préparer la référence* (Prepare reference) → the retained anchor points appear in blue
   on the preview.
4. Adjust if needed (see table), then *▶ Démarrer l'alignement* (Start alignment).

### Settings

| Control (on-screen label) | Default | Purpose |
|---|---|---|
| Local correction (*Correction locale*) | on | off = global shift only |
| Box size (*Boîte*) | 128 px | measurement patch size. Small = tracks fine detail but more fragile; large = robust but smooth |
| Grid step (*Pas de grille*) | 64 px | anchor spacing. Tighter = finer field, longer computation |
| Max local shift (*Décalage local max*) | 12 px | beyond this the point is rejected (outlier guard) |
| Refinement passes (*Passes de raffinement*) | 2 | 2 is enough in most cases; the image is still resampled only once |

The number of retained points is shown live — if it is very low the image is probably too
smooth: reduce the box size and the grid step.

### Output and Siril integration

Aligned frames are written to an **isolated sub-folder** named after the chosen prefix
(default `mp_`), next to the sequence. The isolation is deliberate: Siril's `convert`
command would otherwise mix originals and aligned frames into the same sequence.

Steps chained automatically: *Convertir en séquence Siril* (convert to a Siril sequence),
*Empiler* (stack — `sum`, `med`, `rej 3 3`, `max`, `min`) and *Charger le résultat dans
Siril* (load the result).

---

## ☾ Eclipse mode — Moon (across exposures)

### How it works

Lunar disc detection is **adaptive**, because a corona exposure and a diamond-ring exposure
offer very different cues:

- **Dark disc enclosed** by the corona → centroid of the dark component (accurate, and
  insensitive to corona asymmetry). Confidence 1.00.
- Otherwise (diamond ring, very faint exposures) → **RANSAC circle fit on the limb**,
  refined by Kåsa least squares, robust even on a partial arc. The brightest peak is
  clipped first, otherwise the diamond's glare pulls the detection towards itself.

Frames are then registered by translation so that the centres coincide.

### Steps

1. Exposures from the current Siril sequence load automatically; otherwise
   *⟳ Séquence Siril courante* (current Siril sequence), or *Ajouter* (Add — FITS /
   Sequence / Folder).
2. *↻ Détecter la Lune (toutes)* — Detect the Moon on all frames.
3. **Check the circles**, starting with the low-confidence rows (*Conf.* column: orange
   below 0.15, red below 0.07). This step determines everything downstream.
4. *▶ Aligner sur la Lune* — Align on the Moon.

### Fixing a circle by hand

| Gesture | Effect |
|---|---|
| Click-drag on the preview | moves the centre (committed on release) |
| Shift + wheel | adjusts the radius |
| Wheel | zoom · **Middle-click drag**: pan the view |

Two buttons save a lot of time:

- *Verrouiller au rayon médian* — apply the median detected radius to every frame.
- *⊕ Appliquer le centre affiché à toutes les poses* — the Moon barely moves within a
  bracket: set a perfect centre on one sharp frame, then propagate it to all.

### Output

Files `aligned_00001.fit`, `aligned_00002.fit`… in the exposures' folder (or the folder you
specify), plus the sequence **`aligned_.seq`** if the option is ticked.

Two points that matter downstream:

- Headers are **preserved**, including `EXPTIME` — exactly what FusionHDR needs.
- The output is divided by a **common factor** to fit within [0, 1], so the **relative**
  scale between exposures is preserved — a prerequisite for a radiometrically correct HDR
  merge.

`aligned_.seq` is written directly rather than through `convert`, so that it references
**only** the aligned frames even if other FITS files are present in the folder.

---

## Place in the eclipse workflow

```
ConvertParVitesse  →  SirilJ Align ☾  →  stack (per shutter speed)  →  masters
                          ↓
                   SirilJ Align ☾  (co-register the masters)
                          ↓
                     FusionHDR  →  Corona  →  GHS
```

Eclipse mode is used **twice**:

1. **Before each stack**, on the sub-exposures of a single level. This is essential: if the
   Moon drifts during the stack, its edge is smeared and the master keeps a permanent dark
   bite that no later processing can recover.
2. **On the masters**, to co-register them before the HDR merge.

> When stacking, use `-nonorm`: any normalisation breaks the
> value ∝ exposure × radiance relationship the HDR merge relies on.

See the [full procedure](Procedure_Complete_FR.md) (in French).

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| "Aucune séquence chargée" | no sequence open in Siril | open one, then *⟳ Rafraîchir* |
| "séquence mono-fichier" | FITSEQ or SER | export to individual FITS |
| Few anchor points (Sun) | image too smooth or too noisy | reduce box size and grid step |
| "décalage global aberrant" | frame too different from the reference | change reference, or drop that frame |
| Wrong lunar circle on the diamond ring | detection pulled by the bright point | fix it by click-drag, then *⊕ Appliquer à toutes* |
| Low confidence on several frames | corona too faint to outline the disc | set one good circle and propagate it |

---

## Credits

- Multi-point phase correlation — solar *lucky imaging* principle, in the spirit of AutoStakkert!
- Limb detection — RANSAC circle fit + Kåsa least-squares refinement
- PyQt6 interface, VeraLux conventions
- Licence GPL-3.0-or-later
