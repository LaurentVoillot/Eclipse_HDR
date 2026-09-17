# SPDX-License-Identifier: GPL-3.0-or-later
"""Entrées/sorties images du pipeline autonome.

Format de travail : **TIFF float32, linéaire**, métadonnées en JSON dans le
tag TIFF `ImageDescription`. Photoshop lit le float32 et ignore le tag ; nos
outils y retrouvent le temps de pose, le cercle lunaire et l'historique des
étapes. Les métadonnées voyagent donc *avec* l'image, sans fichier annexe à
perdre.

Convention mémoire : **(H, W, C) float32**, C ∈ {1, 3} — celle de TIFF et
d'OpenCV. Attention : les scripts Siril utilisent la convention FITS
**(C, H, W)** ; `from_fits` / `to_fits` font la transposition.

Toutes les écritures sont **atomiques** (fichier temporaire puis `os.replace`).
Une étape interrompue ne laisse jamais un fichier tronqué qui aurait l'air
valide — c'est ce qui permet de reprendre un traitement sans se demander si
le dernier fichier écrit est complet.

Les intermédiaires sont écrits **non compressés** par défaut : c'est ce qui
rend possible la lecture par bandes (`open_memmap`) dont dépend l'empilement
à mémoire bornée. La compression reste disponible pour les livrables.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any, Optional

import numpy as np
import tifffile

__all__ = [
    "write", "read", "read_meta", "open_memmap",
    "write_mask", "read_mask",
    "export16", "autostretch", "from_fits", "to_fits",
]

# Marqueur de nos propres métadonnées dans ImageDescription.
META_KEY = "eclipse_hdr"
META_VERSION = 1


# ── Écriture / lecture ────────────────────────────────────────────────────────
def _atomic(path: Path, writer) -> Path:
    """Exécute `writer(tmp)` puis renomme sur `path`. Le temporaire est dans le
    même dossier, donc le renommage est atomique (même système de fichiers)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    try:
        writer(tmp)
        os.replace(tmp, path)
    except BaseException:
        tmp.unlink(missing_ok=True)   # y compris sur KeyboardInterrupt
        raise
    return path


def write(path, arr: np.ndarray, meta: Optional[dict] = None, *,
          dtype: str = "float32", compress: bool = False) -> Path:
    """Écrit un TIFF de travail. `meta` est sérialisé en JSON dans le tag
    ImageDescription.

    `dtype="uint16"` stocke les valeurs de `[0,1]` sur seize bits entiers et
    note le facteur dans `scale` ; `read` rend alors les mêmes flottants. À
    réserver aux données qui *sont* déjà quantifiées sur seize bits — la sortie
    du décodage RAW, par exemple, où le float32 ne ferait que rembourrer des
    entiers. Écrire hors de `[0,1]` en entier est refusé plutôt que rogné.

    `compress=True` divise la taille par ~2 mais **interdit la lecture par
    bandes** : à réserver aux fichiers qu'on ne réempilera pas.
    """
    a = np.asarray(arr)
    if a.ndim == 3 and a.shape[2] == 1:
        a = a[:, :, 0]

    payload = {META_KEY: META_VERSION}
    if meta:
        payload.update(meta)

    dt = np.dtype(dtype)
    if dt.kind == "u":
        lo, hi = float(a.min()), float(a.max())
        if lo < -1e-6 or hi > 1.0 + 1e-6:
            raise ValueError(
                f"écriture en {dt.name} demandée mais les valeurs sortent de "
                f"[0,1] ({lo:g} à {hi:g}) — elles seraient rognées ; "
                f"utilisez dtype='float32'")
        full = float(np.iinfo(dt).max)
        payload["scale"] = full
        a = np.rint(np.clip(a, 0.0, 1.0) * full).astype(dt)
    else:
        # `scale` décrit le stockage, pas la donnée. Le laisser survivre à un
        # aller-retour de métadonnées (16 bits → traitement → float32) ferait
        # relire le fichier 65535 fois trop sombre.
        payload.pop("scale", None)
        a = np.ascontiguousarray(a, dtype=dt)

    desc = json.dumps(payload, ensure_ascii=False, separators=(",", ":"))

    photometric = "rgb" if (a.ndim == 3 and a.shape[2] == 3) else "minisblack"
    kw: dict[str, Any] = dict(photometric=photometric, description=desc,
                              metadata=None)   # None : pas de JSON tifffile en plus
    if compress:
        kw["compression"] = "zlib"

    return _atomic(path, lambda p: tifffile.imwrite(str(p), a, **kw))


def read(path) -> tuple[np.ndarray, dict]:
    """→ (données float32 (H,W[,C]), métadonnées).

    Les fichiers écrits en entiers sont ramenés à leur échelle physique : le
    reste du pipeline ne voit que des flottants dans `[0,1]`, quel que soit le
    format de stockage.
    """
    with tifffile.TiffFile(str(path)) as tf:
        arr = tf.asarray()
        meta = _parse_desc(tf.pages[0].description)
    out = np.asarray(arr, dtype=np.float32)
    scale = meta.get("scale")
    if scale:
        out = out / np.float32(scale)
    return out, meta


def read_meta(path) -> dict:
    """Métadonnées seules — ne lit pas les pixels."""
    with tifffile.TiffFile(str(path)) as tf:
        return _parse_desc(tf.pages[0].description)


def _parse_desc(desc) -> dict:
    if not desc:
        return {}
    try:
        d = json.loads(desc)
    except (ValueError, TypeError):
        return {}
    return d if isinstance(d, dict) else {}


def open_memmap(path) -> np.ndarray:
    """Vue mémoire-mappée en lecture seule, pour traiter par bandes sans tout
    charger. Lève `ValueError` si le TIFF est compressé ou non contigu.

    Les valeurs sont **brutes** : pour un fichier entier, il reste à diviser
    par `scale` (voir `open_scaled`)."""
    return tifffile.memmap(str(path), mode="r")


def open_scaled(path) -> tuple[np.ndarray, float]:
    """→ (vue mémoire-mappée, facteur d'échelle à appliquer).

    Le facteur vaut 1 pour un fichier flottant. Il permet de traiter par bandes
    sans se soucier du format de stockage.
    """
    return open_memmap(path), float(read_meta(path).get("scale") or 1.0)


# ── Masques (saturation, rejet…) ──────────────────────────────────────────────
def write_mask(path, mask: np.ndarray) -> Path:
    """Masque en uint8 (0–255 = fraction 0–1), compressé. Un masque de
    saturation est presque vide : il tient en quelques kilo-octets."""
    m = np.asarray(mask)
    u = (m.astype(np.uint8) * 255) if m.dtype == bool else \
        np.clip(np.rint(np.asarray(m, np.float32) * 255.0), 0, 255).astype(np.uint8)
    return _atomic(path, lambda p: tifffile.imwrite(
        str(p), np.ascontiguousarray(u), photometric="minisblack",
        compression="zlib", metadata=None))


def read_mask(path) -> np.ndarray:
    """→ fraction float32 dans [0,1]."""
    return tifffile.imread(str(path)).astype(np.float32) / 255.0


# ── Étirement et export 16 bits ───────────────────────────────────────────────
def _mtf(x, m, lo, hi):
    """Midtone Transfer Function (identique aux scripts Siril)."""
    dist = hi - lo
    if dist < 1e-9:
        return np.where(x > lo, 1.0, 0.0).astype(np.float32)
    xp = np.clip((x - lo) / dist, 0.0, 1.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        y = ((m - 1.0) * xp) / ((2.0 * m - 1.0) * xp - m)
    return np.nan_to_num(y, nan=0.0, posinf=1.0, neginf=0.0).astype(np.float32)


def autostretch(arr: np.ndarray, *, preserve_color: bool = True) -> np.ndarray:
    """Autostretch type Siril (médiane portée à 0,25 ; pied à −2,8 σ robustes).

    `preserve_color` applique la transformation à la **luminance** puis
    réapplique le rapport aux canaux : les couleurs de la couronne restent
    celles de l'image linéaire au lieu d'être désaturées par un étirement
    canal par canal.
    """
    x = np.asarray(arr, np.float32)
    lum = x if x.ndim == 2 else x.mean(axis=2)
    mn, mx = float(lum.min()), float(lum.max())
    ln = (lum - mn) / (mx - mn) if mx > mn else np.zeros_like(lum)

    st = max(1, ln.size // 500_000)
    samp = ln.flatten()[::st]
    med = float(np.median(samp))
    mad = float(np.median(np.abs(samp - med))) * 1.4826 or 0.001
    c0 = max(0.0, med - 2.8 * mad)
    mt = float(_mtf(np.float32(med - c0), 0.25, 0.0, 1.0))
    out_l = np.clip(_mtf(ln, mt, c0, 1.0), 0, 1)

    if x.ndim == 2:
        return out_l
    if not preserve_color:
        xn = (x - mn) / (mx - mn) if mx > mn else np.zeros_like(x)
        return np.clip(_mtf(xn, mt, c0, 1.0), 0, 1)
    gain = out_l / np.maximum(ln, 1e-9)
    xn = (x - mn) / (mx - mn) if mx > mn else np.zeros_like(x)
    return np.clip(xn * gain[..., None], 0, 1)


def _asinh(arr: np.ndarray, k: float = 100.0) -> np.ndarray:
    x = np.asarray(arr, np.float32)
    mx = float(x.max())
    xn = x / mx if mx > 0 else x
    return np.clip(np.arcsinh(k * xn) / np.arcsinh(k), 0, 1).astype(np.float32)


def export16(path, arr: np.ndarray, *, transform: str = "autostretch",
             asinh_k: float = 100.0, scale: Optional[float] = None,
             meta: Optional[dict] = None) -> Path:
    """Exporte en **TIFF 16 bits** pour Photoshop, DxO, Affinity…

    `transform` :
      • `autostretch` (défaut) — étirement automatique, couleurs préservées ;
      • `asinh`      — étirement doux paramétrable, garde le cœur non écrasé ;
      • `linear`     — mise à l'échelle seule (`scale`, sinon 1/max).

    Le défaut n'est **pas** `linear` à dessein : un HDR linéaire couvre cinq
    décades ; ramené à 16 bits sans étirement, la couronne externe tombe sur
    quelques dizaines de niveaux et le fichier s'ouvre en rectangle noir.
    """
    x = np.asarray(arr, np.float32)
    if transform == "autostretch":
        y = autostretch(x)
    elif transform == "asinh":
        y = _asinh(x, asinh_k)
    elif transform == "linear":
        s = scale if scale else (1.0 / max(float(x.max()), 1e-12))
        y = np.clip(x * s, 0, 1)
    else:
        raise ValueError(f"transformation inconnue : {transform!r}")

    u16 = np.rint(np.clip(y, 0, 1) * 65535.0).astype(np.uint16)
    m = dict(meta or {})
    m.pop("scale", None)          # cf. `write` : décrit le stockage d'origine
    m["linear"] = (transform == "linear")
    m["export"] = {"transform": transform,
                   **({"asinh_k": asinh_k} if transform == "asinh" else {}),
                   **({"scale": scale} if transform == "linear" and scale else {})}
    desc = json.dumps({META_KEY: META_VERSION, **m}, ensure_ascii=False)
    photometric = "rgb" if u16.ndim == 3 and u16.shape[2] == 3 else "minisblack"
    return _atomic(path, lambda p: tifffile.imwrite(
        str(p), np.ascontiguousarray(u16), photometric=photometric,
        description=desc, metadata=None, compression="zlib"))


# ── Interopérabilité Siril (FITS) ─────────────────────────────────────────────
def from_fits(path) -> tuple[np.ndarray, dict]:
    """Importe un FITS Siril → ((H,W[,3]) float32, métadonnées).

    Permet d'injecter des masters produits par la chaîne Siril à n'importe
    quelle étape, donc de comparer les deux chaînes sur les mêmes données.
    """
    from astropy.io import fits
    with fits.open(str(path)) as hdul:
        data = np.asarray(hdul[0].data, dtype=np.float32)
        hdr = hdul[0].header
    if data.ndim == 3 and data.shape[0] in (3, 4):
        data = data[:3].transpose(1, 2, 0)
    meta: dict[str, Any] = {}
    for k, key in (("exptime", "EXPTIME"), ("exptime", "EXPOSURE")):
        if key in hdr and "exptime" not in meta:
            meta[k] = float(hdr[key])
    if all(k in hdr for k in ("MOONX", "MOONY", "MOONR")):
        meta["moon"] = {"x": float(hdr["MOONX"]), "y": float(hdr["MOONY"]),
                        "r": float(hdr["MOONR"])}
    for k, key in (("camera", "INSTRUME"), ("stackcnt", "STACKCNT"),
                   ("livetime", "LIVETIME")):
        if key in hdr:
            meta[k] = hdr[key]
    meta["imported_from"] = str(path)
    return data, meta


def to_fits(path, arr: np.ndarray, meta: Optional[dict] = None) -> Path:
    """Exporte en FITS float32 pour reprendre le traitement dans Siril."""
    from astropy.io import fits
    a = np.asarray(arr, np.float32)
    if a.ndim == 3:
        a = a.transpose(2, 0, 1)          # (H,W,C) → (C,H,W), convention FITS
    hdr = fits.Header()
    m = meta or {}
    if m.get("exptime"):
        hdr["EXPTIME"] = (float(m["exptime"]), "s, temps de pose")
        hdr["EXPOSURE"] = (float(m["exptime"]), "s, temps de pose")
    moon = m.get("moon")
    if moon:
        hdr["MOONX"] = (float(moon["x"]), "px, centre Lune")
        hdr["MOONY"] = (float(moon["y"]), "px, centre Lune")
        hdr["MOONR"] = (float(moon["r"]), "px, rayon Lune")
    hdr["HISTORY"] = "eclipse_hdr"
    return _atomic(path, lambda p: fits.PrimaryHDU(
        np.ascontiguousarray(a), header=hdr).writeto(str(p), overwrite=True))
