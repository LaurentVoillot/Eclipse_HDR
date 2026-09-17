# SPDX-License-Identifier: GPL-3.0-or-later
"""Temps de pose : lecture EXIF et regroupement par vitesse.

Porté de `scripts/ConvertParVitesse.py` — même logique de tolérance et mêmes
noms de groupes, pour que les deux chaînes produisent une arborescence
identique et restent comparables.

Le temps de pose est la donnée dont dépend toute la fusion HDR : sans lui, la
relation *valeur ∝ pose × luminance* n'est pas reconstructible. Il est lu une
seule fois ici, puis porté dans les métadonnées de chaque TIFF jusqu'au bout.
"""

from __future__ import annotations

import re
import shutil
import subprocess
from fractions import Fraction
from pathlib import Path
from typing import Optional

__all__ = ["RAW_EXT", "read_exposure", "speed_label", "speed_text",
           "group_by_speed", "scan_folder", "natural_key"]

# Extensions lisibles par libraw, plus les formats matriciels courants.
RAW_EXT = {
    ".cr2", ".cr3", ".crw", ".nef", ".nrw", ".arw", ".srf", ".sr2", ".dng",
    ".orf", ".pef", ".ptx", ".raf", ".rw2", ".raw", ".rwl", ".mrw", ".srw",
    ".x3f", ".3fr", ".fff", ".iiq", ".mos", ".kdc", ".dcr", ".erf",
}


def natural_key(name: str):
    return [int(c) if c.isdigit() else c.lower() for c in re.split(r"(\d+)", str(name))]


# ── Lecture du temps de pose ──────────────────────────────────────────────────
def _exposure_exifread(path: str):
    """→ (num, den, secondes) ou None."""
    try:
        import exifread
    except ImportError:
        return None
    try:
        with open(path, "rb") as fh:
            tags = exifread.process_file(fh, details=False)
    except OSError:
        return None
    for key in ("EXIF ExposureTime", "Image ExposureTime"):
        t = tags.get(key)
        if t is not None and getattr(t, "values", None):
            try:
                r = t.values[0]
                num, den = int(r.num), int(r.den)
                if num > 0 and den > 0:
                    return num, den, num / den
            except (AttributeError, TypeError, ValueError):
                pass
    return None


def _exposure_exiftool(path: str, exe: str):
    """Repli pour les conteneurs non-TIFF (CR3…) → (num, den, secondes)."""
    try:
        out = subprocess.run([exe, "-s3", "-n", "-ExposureTime", str(path)],
                             capture_output=True, text=True, timeout=30)
        v = out.stdout.strip()
        if not v:
            return None
        sec = float(v)
        if sec <= 0:
            return None
        fr = Fraction(sec).limit_denominator(100000)
        return fr.numerator, fr.denominator, sec
    except (OSError, ValueError, subprocess.SubprocessError):
        return None


def read_exposure(path, exiftool: Optional[str] = None):
    """→ (num, den, secondes) ou None. Essaie exifread puis exiftool."""
    r = _exposure_exifread(str(path))
    if r:
        return r
    exe = exiftool or shutil.which("exiftool")
    return _exposure_exiftool(path, exe) if exe else None


# ── Libellés ──────────────────────────────────────────────────────────────────
def speed_label(sec: float, num: Optional[int] = None,
                den: Optional[int] = None) -> str:
    """Nom de dossier : 1/160 s → « 1-160 », 2 s → « 2s », 1,3 s → « 1s3 ».
    Jamais de « / » ni de « . » — noms portables sur tous les systèmes."""
    if sec >= 1.0:
        if abs(sec - round(sec)) < 1e-6:
            return f"{int(round(sec))}s"
        txt = f"{sec:.2f}".rstrip("0").rstrip(".")
        return txt.replace(".", "s")
    if num == 1 and den and den > 1:
        return f"1-{den}"
    inv = 1.0 / sec
    return f"1-{int(round(inv))}" if inv >= 1 else f"1-{inv:.3g}".replace(".", "s")


def speed_text(sec: float, num: Optional[int] = None,
               den: Optional[int] = None) -> str:
    """Libellé lisible : « 1/160 s », « 2 s »."""
    if sec >= 1.0:
        return f"{sec:g} s"
    if num == 1 and den:
        return f"1/{den} s"
    return f"1/{round(1.0 / sec)} s"


# ── Regroupement ──────────────────────────────────────────────────────────────
def group_by_speed(entries, tol_pct: float = 3.0):
    """Regroupe des poses par temps d'exposition.

    Deux poses sont dans le même groupe si leur écart relatif est sous la
    tolérance (3 % par défaut : fusionne les arrondis APEX du type 1/1000 vs
    1/1024, sans mélanger deux crans de 1/3 de diaph, distants de 26 %).

    → ([{label, sec, num, den, files:[…]}], inconnus)
    """
    known = [e for e in entries if e.get("sec")]
    unknown = [e for e in entries if not e.get("sec")]
    known.sort(key=lambda e: (e["sec"], natural_key(e["path"])))
    groups = []
    for e in known:
        if groups and abs(e["sec"] - groups[-1]["sec_ref"]) <= \
                (tol_pct / 100.0) * groups[-1]["sec_ref"]:
            groups[-1]["files"].append(e)
            continue
        groups.append({"sec_ref": e["sec"], "files": [e]})

    out = []
    for g in groups:
        mid = g["files"][len(g["files"]) // 2]      # pose médiane du groupe
        out.append({
            "label": speed_label(mid["sec"], mid.get("num"), mid.get("den")),
            "sec": mid["sec"], "num": mid.get("num"), "den": mid.get("den"),
            "files": g["files"],
        })
    return out, unknown


def scan_folder(raw_dir, tol_pct: float = 3.0, exiftool: Optional[str] = None,
                progress=None):
    """Inventorie un dossier de RAW → (groupes, inconnus).

    Ne convertit rien : lit seulement les EXIF. `progress(i, n, nom)` est
    appelé à chaque fichier, pour l'affichage.
    """
    raw_dir = Path(raw_dir)
    files = sorted((f for f in raw_dir.iterdir()
                    if f.is_file() and f.suffix.lower() in RAW_EXT),
                   key=lambda f: natural_key(f.name))
    exe = exiftool or shutil.which("exiftool")
    entries = []
    for i, f in enumerate(files):
        if progress:
            progress(i, len(files), f.name)
        r = read_exposure(f, exe)
        e = {"path": str(f), "name": f.name}
        if r:
            e["num"], e["den"], e["sec"] = r
        entries.append(e)
    return group_by_speed(entries, tol_pct)
