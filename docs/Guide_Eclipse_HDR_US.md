# Guide — Total Solar Eclipse HDR Processing

**From RAW to a finished corona in Siril.**
Full chain: capture → conversion → per-exposure alignment/stacking → HDR merge → coronal enhancement → finishing.

> **Note:** the `FusionHDR` and `Corona` scripts have a **French** interface. This guide is in English; their French button labels are shown in **bold** with an English gloss.

---

## ⚠ Safety — eyes & gear (read first)

**Outside totality, the Sun is dangerous to both your eyes and your sensor.** The photosphere is visible during every partial phase (before C2, after C3): looking at it — or pointing an optic at it — without protection causes **irreversible, painless retinal burns**, and can destroy the sensor or start a fire.

| Phase | Eyes | Optics / camera |
|---|---|---|
| **Partial phases** (before C2, after C3) | **Certified ISO 12312-2 eclipse glasses** only. Never with the naked eye, nor through a viewfinder / binoculars / telescope without a filter. | **Certified solar filter** (e.g. Baader AstroSolar) **in front** of the lens, at all times. |
| **Totality** (between C2 and C3, **only**) | Naked-eye viewing **allowed** — the only safe moment. | Shooting **without a filter** allowed. |

- **At C2 / end of the diamond ring**: remove the filter **only** once Baily's beads and the diamond ring are gone and totality is fully established.
- **At C3 / reappearance of the diamond ring**: the photosphere returns **suddenly**. **Put the filter back and stop all direct viewing the instant the diamond ring reappears** — do not wait: a fraction of a second can injure.
- Know your precise **contact times (C2, C3)**, keep the filter **ready to refit**, and warn everyone present. When in doubt: **filter on**.

> Eclipse glasses protect **your eyes**, not the camera; the lens filter protects **the camera**, not your eyes if you look elsewhere. Both are required during the partial phases.

---

## Tools in the chain

| Step | Tool |
|---|---|
| RAW → FITS conversion | **Siril** (*Conversion* tab) |
| Per-level alignment + stacking | **Siril** (*Registration* + *Stacking*) or **SolarAlign** |
| Co-aligning the masters | **Siril** (DFT registration) |
| 32-bit HDR merge | **FusionHDR** |
| Coronal reveal | **Corona** |
| Curves / finishing | **Siril** (Histogram / Asinh / GHS) |

> Software prerequisites: Siril 1.3+, and the `FusionHDR.py` / `Corona.py` scripts in Siril's scripts folder.

---

## Phase 0 — Capture (eclipse day)

The quality of the result is decided here.

- **Equatorial mount** tracking, well polar-aligned → the Sun stays framed, alignment is nearly pure translation (no field rotation).
- **Manual focus** locked (set it beforehand on the lunar limb or a bright star), then **tape the focus ring**. Don't touch zoom/focus afterwards.
- **Manual exposure**: fixed ISO (100–400), fixed aperture, vary only the **shutter speed**.
- **Wide bracketing**: from ~**1/4000 s to ~2–4 s**, in ~1-stop steps → **~12–15 levels**. The corona spans a huge brightness ratio; you need the whole range.
- **RAW** mandatory (linear data). Continuous bursts, intervalometer or in-camera bracketing, **identical framing**.
- **Solar filter**: remove it **at C2** (totality fully established, diamond ring gone) and **refit it at C3 the moment the diamond ring reappears** — rehearse this gesture before the day. See the **⚠ Safety** section: it is vital for your eyes **and** the sensor.
- Plan the specific moments: **Baily's beads / diamond ring** (very short exposures) and **prominences** (short exposures), just after C2 and before C3.
- Anti-vibration: self-timer, electronic shutter or mirror lock-up.

> ⚠ **The Moon moves across the corona** (~0.5″/s, i.e. ~2′ over 4 min). See "Pitfalls": this drives the alignment choice.

---

## Phase 1 — Convert RAW to FITS (Siril)

1. **Conversion** tab.
2. Drop all RAW files (or by group — see Phase 2).
3. Output **FITS**, tick **Debayer** for a **color** camera (OSC/DSLR).
4. **Convert**. Siril produces linear FITS (`EXPTIME`/EXIF is preserved — useful for FusionHDR).

> Keep the data **linear**; don't stretch anything at this stage.

---

## Phase 2 — Sort by exposure level

The HDR merge combines **one master per brightness**, so group the frames by shutter speed.

- Sort the RAW files into **per-exposure folders**: `1_2000/`, `1_500/`, `1_125/`, … `2s/` (from the EXIF shutter speed).
- Convert **each folder** into its own FITS sequence (Phase 1 per folder).

Result: one FITS sequence per exposure level.

---

## Phase 3 — Align + stack each level

For **each** sequence (each brightness):

**Option A — Siril native (recommended)**
1. **Registration** tab → method **"Image Pattern Alignment (DFT)"** (phase correlation, translation). It locks onto the dominant high-contrast pattern: **lunar disk + inner corona**. ("Global Star Alignment" won't work: no stars.)
2. **Register**.
3. **Stacking** tab → **Average with rejection** (Winsorized/Sigma, e.g. 3/3) to remove planes, satellites, cloud gusts.
4. **Stack** → one **master** for this level (e.g. `master_1_500.fit`).

**Option B — SolarAlign script**
Load the level's sequence in Siril, run **SolarAlign** (StackReg *Corps rigide* / Rigid body, or ECC *Euclidien* / Euclidean), tick *Convertir* (Convert) + *Empiler* (Stack) → it aligns and stacks in one go. Handy if native registration struggles.

> **Long vs short exposures tip:**
> - **Long exposures** (faint outer corona) → stack many (reduces noise); the overexposed lunar edge will be masked anyway.
> - **Short exposures** (prominences, inner corona) → use a **short time window** or single best frames (prominences evolve, the Moon drifts).

---

## Phase 4 — Co-align the masters

The masters must **overlay** before the HDR merge (otherwise edge artifacts).

1. Put all masters in one folder, convert them into **a single sequence**.
2. **Registration** → **DFT** (the lunar disk, common to every exposure, anchors the alignment despite the different brightnesses).
3. Export/save the **registered** frames (prefix `r_`).

> If the mount tracked perfectly with strictly identical framing, the masters may already be aligned; check, otherwise do this phase.

---

## Phase 5 — HDR merge (FusionHDR)

> **No bracketing?** This phase is **optional**. If you have a single exposure (one stacked master), **skip Phase 5** and go straight to **Corona** (Phase 6).

1. **Scripts → FusionHDR**.
2. **Ajouter…** (Add) the registered masters (one per exposure).
3. Check the **Expo (s)** column (read from the header). If missing/identical → **Estimer les expositions (recouvrement)** (Estimate exposures from overlap).
4. **Mode**: *Remplacement par seuil* (Threshold replacement — default, most predictable).
5. **Pleine échelle** (Full scale): *Auto*. **Seuil saturation** (Saturation threshold) ~0.95.
6. **Soustraire le fond** (Subtract background): **unchecked** for eclipse (sky ≈ black).
7. **Normaliser la sortie** (Normalize output): checked.
8. **Fusionner** (Merge) → preview. **Enregistrer + charger dans Siril** (Save + load) → a **32-bit linear FITS** HDR.

---

## Phase 6 — Reveal the corona (Corona)

On the HDR loaded in Siril:

1. **Scripts → Corona**.
2. **Center**: click **Auto**, then fine-tune by **clicking** the disk center. Set the **lunar radius** so the **blue circle** matches the Moon's edge.
3. **Retrait du gradient radial** (Radial gradient removal): **enabled** (the key step — it flattens the falloff and brings out the streamers). Smoothing ~8.
4. **Appliquer** (Apply) → the corona already appears structured.
5. Optional **Larson–Sekanina**: Δ angle ~3–8°, Δ radial 0, amount ~1 → crisps the **streamers** (radial structures). Increase carefully (amplifies noise).
6. Optional **Unsharp**: sigma 2–4, amount ~1 → fine detail.
7. **Noircir le disque lunaire** (Blacken the lunar disk): checked.
8. **Enregistrer + charger dans Siril** (Save + load).

> Iterate: adjust center / radius / parameters and **Apply** again until the look is right.

---

## Phase 7 — Finishing in Siril (the curves)

Corona's output is flattened: now for the aesthetics.

- **Histogram Transformation** (stretch) or **Asinh** to reveal the corona's soft dynamic range.
- **Generalized Hyperbolic Stretch (GHS)** for fine local-contrast control without blowing out the center.
- **Color**: white balance / *Photometric Color Calibration* if relevant; moderate saturation (the corona has subtle color structure).
- Light **noise reduction** on the faint regions if needed.
- Final export: **16-bit TIFF** or **PNG**.

---

## Pitfalls & variants

- **Moon vs corona**: the corona is fixed (to the Sun), the Moon drifts. For a long totality, process **two versions** — one **corona-aligned** (long exposures, outer corona) and one **Moon-aligned** (short exposures, prominences/inner corona) — then **composite** (the corona + the Moon/prominences captured at the right instant). For a first pass or a short totality, DFT alignment on the lunar disk is enough.
- **Prominences & chromosphere**: come from **short** exposures at a precise instant (just after C2 / before C3), **not stackable** across the whole totality. Add them as a layer/composite.
- **Lunar disk**: often replaced by a clean black disk, or by an **earthshine Moon** shot separately (long exposure).
- **Color**: FusionHDR and Corona support color (radial model from luminance, applied per channel); keep the white balance for the end of the chain.
- **No HDR (single exposure)**: the chain also works **without bracketing**. Stack your frames (Phase 3), then apply **Corona** directly to the master, **skipping Phase 5** (FusionHDR). Corona reveals the corona on **any** image — HDR or not — as long as the center isn't heavily saturated. (FusionHDR itself needs at least two exposures of different brightness.)
- **Memory**: FusionHDR merges **one image at a time** (scales to large files); Corona loads the whole image.

---

## Quick recap

```
RAW (bracketing, RAW, tracking mount, filter removed at C2)
   │  Siril: Conversion (+ Debayer if color) → linear FITS
   ▼
Sort by exposure (folders)
   │  per level: DFT Registration + Stacking (rejection) → 1 master/exposure
   ▼  (or SolarAlign: align + stack)
Co-align the masters (Siril DFT) → registered masters
   │
   ▼  FusionHDR: threshold replacement, auto scale → 32-bit HDR   (optional — skip if no bracketing)
   │
   ▼  Corona: center + radial gradient removal (+ LS + unsharp) → corona revealed
   │
   ▼  Siril: Histogram / Asinh / GHS, color, export TIFF/PNG
RESULT
```

---

*This same FusionHDR + Corona chain is reusable beyond the eclipse: FusionHDR for lunar (earthshine, lunar eclipse), planetary (planet + moons) and deep-sky (the bright core of M42…); Corona stays specific to the solar corona.*
