##############################################
# EclipseComposite — couches expérimentales
# Lune (clair de Terre) · Protubérances · Étoiles
# Version 0.1.0 — outil exploratoire
##############################################
#
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Composition ASSUMÉE comme esthétique (plus radiométrique) : on injecte
# dans le HDR des couches issues des poses alignées, chacune avec son
# niveau relatif. C'est la méthode des grands composites d'éclipse
# (lumière cendrée dans le disque, protubérances colorées, étoiles).
#
# Crédits
# -------
#   • Détection du disque : Otsu + composante sombre enclose (cf. Corona)
#   • Interface : PyQt6 / conventions VeraLux

"""
EclipseComposite — ajoute au HDR fusionné trois couches optionnelles :

☾ Lune (clair de Terre)
    Empilement médian des N poses les PLUS LONGUES restreint au disque
    lunaire (les poses sont alignées sur la Lune → empilement net), avec
    retrait du halo interne (diffusion de la couronne qui déborde sur le
    bord du disque). Niveau réglable en % de la couronne interne.

🔥 Protubérances
    Anneau au limbe pris des N poses les PLUS COURTES (non saturées là),
    couleur préservée. Évite les protubérances cramées en blanc.

✶ Étoiles
    Détection de sources ponctuelles (DoG + seuil kσ + forme PSF) sur
    l'empilement des poses longues, HORS de la couronne, avec option de
    confirmation croisée (présence dans deux demi-empilements) pour
    éliminer les pics de bruit. Rehaussement local uniquement.

Entrées : le HDR (image courante de Siril ou fichier) + les poses
alignées (celles données à FusionHDR, avec EXPTIME dans les en-têtes).
Sortie : FITS 32 bits linéaire, à finir au GHS dans Siril.
"""

import sys
import os
import re
import tempfile

try:
    import sirilpy as s
except ImportError:
    print("Erreur : module sirilpy introuvable.")
    sys.exit(1)

s.ensure_installed("numpy", "astropy", "scipy", "PyQt6")

import numpy as np
from pathlib import Path
from typing import Optional

from astropy.io import fits
from scipy.ndimage import (gaussian_filter, gaussian_filter1d, label,
                           binary_closing, maximum_filter)

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
    QLabel, QPushButton, QCheckBox, QSpinBox, QDoubleSpinBox, QProgressBar,
    QTextEdit, QFileDialog, QGroupBox, QGraphicsView, QGraphicsScene,
    QTableWidget, QTableWidgetItem, QHeaderView, QAbstractItemView
)
from PyQt6.QtCore import Qt, QThread, pyqtSignal, QTimer
from PyQt6.QtGui import QImage, QPixmap, QPainter, QColor, QPen

VERSION = "0.1.0"

DARK_SS = """
QWidget            { background:#2b2b2b; color:#d4d4d4; font-size:11px; }
QMainWindow        { background:#2b2b2b; }
QPushButton        { background:#3c3c3c; color:#d4d4d4; border:1px solid #555;
                     border-radius:4px; padding:4px 10px; }
QPushButton:hover  { background:#4a4a4a; }
QPushButton:disabled { color:#666; border-color:#444; }
QPushButton#BtnGo  { background:#1e6b1e; color:#fff; font-weight:bold;
                     border:1px solid #2a9b2a; padding:6px 12px; }
QPushButton#BtnGo:hover    { background:#248c24; }
QPushButton#BtnGo:disabled { background:#333; color:#666; border-color:#444; }
QSpinBox, QDoubleSpinBox { background:#1e1e1e; border:1px solid #555;
                            border-radius:3px; padding:2px 4px; color:#d4d4d4; }
QCheckBox::indicator { width:14px; height:14px; border:1px solid #666;
                        border-radius:2px; background:#1e1e1e; }
QCheckBox::indicator:checked { background:#4caf50; border-color:#4caf50; }
QGroupBox          { border:1px solid #555; border-radius:4px;
                     margin-top:8px; padding-top:4px; }
QGroupBox::title   { subcontrol-origin:margin; subcontrol-position:top left;
                     padding:0 4px; color:#aaa; }
QProgressBar       { border:1px solid #555; border-radius:3px; background:#1e1e1e;
                     text-align:center; color:#d4d4d4; }
QProgressBar::chunk { background:#2a6f2a; }
QTextEdit          { background:#1a1a1a; border:1px solid #444; color:#b4b4b4;
                     font-family:monospace; font-size:10px; }
QGraphicsView      { background:#111; border:1px solid #444; }
QTableWidget       { background:#1e1e1e; gridline-color:#3a3a3a; color:#d4d4d4; }
QTableWidget::item:selected { background:#3a5f3a; }
QHeaderView::section { background:#333; color:#aaa; padding:3px; border:1px solid #444; }
QScrollBar:vertical   { background:#2b2b2b; width:10px; }
QScrollBar::handle:vertical { background:#555; border-radius:5px; }
"""

FITS_EXT = {".fit", ".fits", ".fts"}


# ── FITS / affichage ──────────────────────────────────────────────────────────
def load_fits(path):
    """FITS → (data float32 valeurs d'origine, header). Couleur → (C,H,W)."""
    with fits.open(path) as hdul:
        hdr = hdul[0].header
        data = hdul[0].data.astype(np.float32)
    if data.ndim == 3 and data.shape[2] in (3, 4) and data.shape[0] not in (3, 4):
        data = data.transpose(2, 0, 1)
    if data.ndim == 3 and data.shape[0] == 4:
        data = data[:3]
    return data, hdr

def get_exptime(hdr):
    for k in ("EXPTIME", "EXPOSURE", "EXPO"):
        if k in hdr:
            try:
                v = float(hdr[k])
                if v > 0:
                    return v
            except Exception:
                pass
    return 0.0

def to_lum(d):
    return d if d.ndim == 2 else d.mean(axis=0)

def save_fits(data, path):
    fits.PrimaryHDU(np.asarray(data, dtype=np.float32)).writeto(path, overwrite=True)

def _mtf(x, m, lo, hi):
    dist = hi - lo
    if dist < 1e-9:
        return np.where(x > lo, 1.0, 0.0).astype(np.float32)
    xp = np.clip((x - lo) / dist, 0.0, 1.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        y = ((m - 1.0) * xp) / ((2.0 * m - 1.0) * xp - m)
    return np.nan_to_num(y, nan=0.0, posinf=1.0, neginf=0.0).astype(np.float32)

def autostretch_display(data):
    chans = [data] if data.ndim == 2 else [data[i] for i in range(data.shape[0])]
    mx = max(float(c.max()) for c in chans) or 1.0
    chans = [np.clip(c / mx, 0, 1) for c in chans]
    lum = chans[0] if len(chans) == 1 else np.mean(chans, axis=0)
    st = max(1, lum.size // 500_000)
    samp = lum.flatten()[::st]
    med = float(np.median(samp))
    mad = float(np.median(np.abs(samp - med))) * 1.4826 or 0.001
    c0 = max(0.0, med - 2.8 * mad)
    mt = float(_mtf(np.float32(med - c0), 0.25, 0.0, 1.0))
    out = np.stack([_mtf(c, mt, c0, 1.0) for c in chans])
    return out[0] if len(chans) == 1 else out

def to_pixmap(disp):
    if disp.ndim == 2:
        g = np.ascontiguousarray(np.flipud((np.clip(disp, 0, 1) * 255).astype(np.uint8)))
        h, w = g.shape
        return QPixmap.fromImage(QImage(g.tobytes(), w, h, w, QImage.Format.Format_Grayscale8))
    rgb = np.flipud((np.clip(disp, 0, 1) * 255).astype(np.uint8).transpose(1, 2, 0))
    rgb = np.ascontiguousarray(rgb[:, :, :3])
    h, w, _ = rgb.shape
    return QPixmap.fromImage(QImage(rgb.tobytes(), w, h, 3 * w, QImage.Format.Format_RGB888))

def smoothstep(lo, hi, x):
    t = np.clip((x - lo) / max(hi - lo, 1e-9), 0.0, 1.0)
    return t * t * (3.0 - 2.0 * t)


# ── Détection du disque lunaire — détecteur ADAPTATIF (repris de SirilJ_Align)
#    Disque sombre enclos par la couronne → centroïde (le plus précis).
#    Sinon (couronne faible dans une direction, le disque « fuit » vers le
#    fond) → ajustement du limbe par RANSAC + Kåsa, qui n'exige qu'un ARC.
#    L'ancien détecteur n'avait que le premier cas et abandonnait souvent.
def _otsu(v):
    vmin, vmax = float(v.min()), float(v.max())
    if vmax <= vmin:
        return vmin
    hist, edges = np.histogram(v.ravel(), bins=256, range=(vmin, vmax))
    hist = hist.astype(np.float64); tot = hist.sum()
    if tot == 0:
        return (vmin + vmax) / 2
    p = hist / tot; c = (edges[:-1] + edges[1:]) / 2
    om = np.cumsum(p); mu = np.cumsum(p * c); mt = mu[-1]
    den = om * (1 - om)
    with np.errstate(divide="ignore", invalid="ignore"):
        sb = np.where(den > 1e-12, (mt * om - mu) ** 2 / den, 0.0)
    return float(c[int(np.argmax(sb))])

def _circumcircle(p1, p2, p3):
    ax, ay = p1; bx, by = p2; cx, cy = p3
    d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    if abs(d) < 1e-6:
        return None
    ux = ((ax*ax+ay*ay)*(by-cy) + (bx*bx+by*by)*(cy-ay) + (cx*cx+cy*cy)*(ay-by)) / d
    uy = ((ax*ax+ay*ay)*(cx-bx) + (bx*bx+by*by)*(ax-cx) + (cx*cx+cy*cy)*(bx-ax)) / d
    return ux, uy, float(np.hypot(ux - ax, uy - ay))

def _fit_circle_kasa(pts):
    x = pts[:, 0]; y = pts[:, 1]
    A = np.column_stack([x, y, np.ones(len(x))])
    b = x * x + y * y
    sol, *_ = np.linalg.lstsq(A, b, rcond=None)
    cx = sol[0] / 2.0; cy = sol[1] / 2.0
    return cx, cy, float(np.sqrt(max(sol[2] + cx * cx + cy * cy, 0.0)))

def _detect_disk_enclosed(a):
    """Centroïde du disque sombre **enclos** par la couronne (précis pour les
    poses lumineuses/moyennes, insensible au halo surexposé). → (cx,cy,r) ou None."""
    bright = a >= _otsu(a)
    try:
        bright = binary_closing(bright, iterations=2)
    except Exception:
        pass
    lbl, n = label(~bright)
    if n == 0:
        return None
    border = set(np.unique(np.concatenate([lbl[0, :], lbl[-1, :], lbl[:, 0], lbl[:, -1]])))
    sizes = np.bincount(lbl.ravel())
    best, bs = 0, 0
    for L in range(1, n + 1):
        if L in border:
            continue
        if sizes[L] > bs:
            best, bs = L, sizes[L]
    if best == 0 or bs < 0.005 * a.size:
        return None
    ys, xs = np.where(lbl == best)
    return float(xs.mean()), float(ys.mean()), float(np.sqrt(bs / np.pi))

def detect_moon(lum, seed=42):
    """Détecte la Lune → (cx, cy, r, confiance) en coords tableau, ou None.

    Adaptatif : disque sombre enclos par la couronne → centroïde (confiance
    1.0) ; sinon (diamant / poses faibles) → limbe par RANSAC (robuste même
    sur un arc partiel)."""
    a = gaussian_filter(lum.astype(np.float32), 2.0)
    cap = float(np.percentile(a, 99.5))     # écrête le diamant (déborde dans le disque)
    if cap > 0:
        a = np.minimum(a, cap)

    disk = _detect_disk_enclosed(a)
    if disk is not None:
        return disk[0], disk[1], disk[2], 1.0

    gy, gx = np.gradient(a)
    gmag = np.hypot(gx, gy)
    H, W = a.shape
    thr = np.percentile(gmag, 99.0)
    ys, xs = np.where(gmag >= thr)
    if len(xs) < 30:
        return None
    pts = np.column_stack([xs, ys]).astype(np.float64)
    rng = np.random.default_rng(seed)
    if len(pts) > 3000:
        pts = pts[rng.choice(len(pts), 3000, replace=False)]
    rmin, rmax = 0.04 * min(H, W), 0.49 * min(H, W)
    best_n, best = 0, None
    for _ in range(1500):
        i = rng.choice(len(pts), 3, replace=False)
        c = _circumcircle(pts[i[0]], pts[i[1]], pts[i[2]])
        if not c:
            continue
        cx, cy, r = c
        if not (rmin <= r <= rmax):
            continue
        d = np.abs(np.hypot(pts[:, 0] - cx, pts[:, 1] - cy) - r)
        nn = int(np.count_nonzero(d < 2.0))
        if nn > best_n:
            best_n, best = nn, (cx, cy, r)
    if best is None:
        return None
    cx, cy, r = best
    for _ in range(2):
        d = np.abs(np.hypot(pts[:, 0] - cx, pts[:, 1] - cy) - r)
        inl = pts[d < 3.0]
        if len(inl) < 10:
            break
        cx, cy, r = _fit_circle_kasa(inl)
    return float(cx), float(cy), float(r), float(best_n / len(pts))


def detect_disk(lum):
    """→ (cx, cy, r) ou None. Conserve la signature d'origine."""
    res = detect_moon(lum)
    return None if res is None else (res[0], res[1], res[2])


# ══════════════════════════════════════════════════════════════════════════════
#  Cœur des couches (fonctions pures, testables hors Qt)
# ══════════════════════════════════════════════════════════════════════════════
def median_stack_radiance(paths_exps):
    """Empilement MÉDIAN des radiances (valeur/EXPTIME) d'une liste de poses.
    Retourne (stack, shape) — mono (H,W) ou couleur (C,H,W)."""
    rads = []
    shape = None
    for p, e in paths_exps:
        d, _ = load_fits(p)
        if shape is None:
            shape = d.shape
        if d.shape != shape:
            raise ValueError(f"dimensions différentes : {Path(p).name}")
        rads.append(d / max(e, 1e-9))
    return np.median(np.stack(rads), axis=0).astype(np.float32)


def remove_inner_halo(layer, cx, cy, R):
    """Retire, À L'INTÉRIEUR du disque, la composante radiale montante vers le
    limbe (diffusion de la couronne qui déborde). Le piédestal de lumière
    cendrée (plat) et les mers (albédo) sont conservés. Par canal."""
    chans = [layer] if layer.ndim == 2 else [layer[i] for i in range(layer.shape[0])]
    H, W = chans[0].shape
    yy, xx = np.indices((H, W)).astype(np.float32)
    r = np.hypot(xx - cx, yy - cy)
    inside = r < R
    nb = max(int(R) + 1, 2)
    ri = np.clip(r.astype(np.int32), 0, nb - 1)
    idx = ri[inside].ravel()
    out = []
    for c in chans:
        vals = c[inside].ravel().astype(np.float64)
        cnt = np.bincount(idx, minlength=nb).astype(np.float64)
        ssum = np.bincount(idx, weights=vals, minlength=nb)
        good = cnt > 0
        prof = np.zeros(nb)
        prof[good] = ssum[good] / cnt[good]
        if good.sum() >= 2:
            ai = np.arange(nb)
            prof = np.interp(ai, ai[good], prof[good])
        prof = gaussian_filter1d(prof, 3.0)
        inner = float(np.median(prof[: max(int(0.5 * R), 1)]))   # piédestal cendré
        halo = np.maximum(prof - inner, 0.0)                     # excès montant au limbe
        corr = halo[ri].astype(np.float32)
        cc = c.copy()
        cc[inside] = np.maximum(c[inside] - corr[inside], 0.0)
        out.append(cc)
    return out[0] if layer.ndim == 2 else np.stack(out)


def blend_moon(base, layer, cx, cy, R, level_frac):
    """Injecte la couche disque dans base, bord doux à l'intérieur du limbe.
    level_frac : médiane du disque amenée à level_frac × P99,9(base)."""
    H, W = (base.shape[-2], base.shape[-1])
    yy, xx = np.indices((H, W)).astype(np.float32)
    r = np.hypot(xx - cx, yy - cy)
    w = 1.0 - smoothstep(max(R - 12.0, 0.0), max(R - 3.0, 1.0), r)   # 1 au centre → 0 au limbe
    lum_l = to_lum(layer)
    inside = r < 0.8 * R
    ref_l = float(np.percentile(lum_l[inside], 90.0)) if inside.any() else 0.0
    ref_b = float(np.percentile(to_lum(base), 99.9))
    gain = (level_frac * ref_b / ref_l) if ref_l > 1e-12 else 0.0
    lay = layer * gain
    if base.ndim == 3 and lay.ndim == 2:
        lay = np.stack([lay] * base.shape[0])
    if base.ndim == 3:
        w = w[None, :, :]
    return (base * (1.0 - w) + lay * w).astype(np.float32), gain


def halpha_weight(rgb, ring, softness=1.0):
    """Poids « rougeur Hα » : 1 sur les protubérances, 0 sur la couronne.

    La couronne K est BLANCHE (diffusion Thomson de la photosphère) alors que
    les protubérances émettent en Hα (656 nm) : leur excès de rouge est le
    discriminant le plus sûr, bien meilleur que l'intensité — c'est pourquoi
    elles se voient sur le disque sombre mais se noient dans la couronne
    brillante, où leur contraste d'intensité est négligeable.
    Retourne un poids [0,1] de même forme que la luminance."""
    if rgb.ndim != 3 or rgb.shape[0] < 3:
        return np.ones(rgb.shape[-2:], np.float32)
    r_, g_, b_ = rgb[0], rgb[1], rgb[2]
    lum = np.maximum((r_ + g_ + b_) / 3.0, 1e-12)
    exc = (r_ - np.maximum(g_, b_)) / lum          # excès de rouge relatif
    v = exc[ring]
    if v.size < 50:
        return np.ones(rgb.shape[-2:], np.float32)
    # seuil robuste : au-dessus de la dispersion du fond (couronne neutre)
    med = float(np.median(v))
    sig = 1.4826 * float(np.median(np.abs(v - med))) or 1e-6
    t = (exc - (med + 1.5 * sig)) / max(2.0 * sig * float(softness), 1e-9)
    t = np.clip(t, 0.0, 1.0)
    return (t * t * (3.0 - 2.0 * t)).astype(np.float32)


def neutralize_corona(img, cx, cy, R, amount, balance=True):
    """Désature la COURONNE en laissant les protubérances rouges intactes.

    La couronne K est de la lumière photosphérique diffusée par les électrons :
    elle est BLANCHE. Toute teinte qu'elle porte est un artefact (balance des
    blancs, dématriçage, étirement). La neutraliser maximise le RAPPORT de
    saturation entre couronne et protubérances — c'est lui qui compte, la
    saturation absolue se remontant ensuite d'un curseur.

    amount : 0 = inchangé, 1 = couronne parfaitement neutre.
    balance : égalise d'abord les canaux sur la couronne (retire la dominante).
    """
    # balance et désaturation sont deux réglages indépendants
    if img.ndim != 3 or img.shape[0] < 3 or (amount <= 1e-6 and not balance):
        return img.astype(np.float32), (1.0, 1.0, 1.0)
    H, W = img.shape[-2], img.shape[-1]
    yy, xx = np.indices((H, W)).astype(np.float32)
    r = np.hypot(xx - cx, yy - cy)
    zone = (r > R + 6.0)                       # couronne, hors disque lunaire
    out = img.astype(np.float32).copy()
    gains = (1.0, 1.0, 1.0)
    if balance and zone.any():
        # Balance des blancs mesurée sur la couronne SEULE, en excluant les
        # pixels rouges (protubérances) pour ne pas les blanchir.
        hw0 = halpha_weight(out, zone)
        neutral = zone & (hw0 < 0.2)
        ref = neutral if int(neutral.sum()) > 500 else zone
        meds = [float(np.median(out[c][ref])) for c in range(3)]
        tgt = float(np.mean(meds))
        if all(m > 1e-9 for m in meds):
            gains = tuple(tgt / m for m in meds)
            for c in range(3):
                out[c] = out[c] * gains[c]
    # Désaturation pondérée : pleine sur le neutre, nulle sur le rouge Hα.
    hw = halpha_weight(out, zone)
    s = (1.0 - float(amount) * (1.0 - hw)).astype(np.float32)
    lum = out.mean(axis=0)
    out = np.maximum(lum[None, :, :] + (out - lum[None, :, :]) * s[None, :, :], 0.0)
    return out.astype(np.float32), gains


def purify_red(img, amount):
    """Rouge pur : retire le VERT des pixels chauds (jaune/orange → rouge).

    Les protubérances émettent en Hα (656 nm), donc essentiellement dans le
    canal R. Le vert qu'on y voit est une fuite : débordement dans les photo-
    sites verts et étalement au dématriçage. Le résultat vire au jaune/orange.

    Critère : la dominance du rouge sur le BLEU — et non sur max(V,B), car
    pour un pixel jaune R ≈ V et un tel test ne le détecterait pas. Là où le
    rouge domine, le vert est ramené vers le bleu : la teinte devient un rouge
    franc, le blanc et le bleu ne bougent pas.

    amount : 0 = inchangé, 1 = vert entièrement ramené au niveau du bleu.
    """
    if img.ndim != 3 or img.shape[0] < 3 or amount <= 1e-6:
        return img.astype(np.float32)
    out = img.astype(np.float32).copy()
    r_, g_, b_ = out[0], out[1], out[2]
    warm = np.clip((r_ - b_) / np.maximum(r_, 1e-12), 0.0, 1.0)
    w = (warm * warm * (3.0 - 2.0 * warm)).astype(np.float32)   # transition douce
    excess = np.maximum(g_ - b_, 0.0)            # vert au-dessus du neutre
    out[1] = g_ - float(amount) * w * excess
    return np.maximum(out, 0.0).astype(np.float32)


def blend_prominences(base, short_rad, cx, cy, R, thickness, level_frac,
                      chroma=1.0, ha_sel=0.0, purify=0.0):
    """AJOUTE uniquement les protubérances de la pose courte (couleur préservée).

    Dans l'anneau, la pose courte contient partout la couronne interne : la
    remplacer en bloc créerait un anneau blanc uniforme. On extrait donc le
    DÉTAIL au-dessus de la ligne de base azimutale (médiane par rayon, par
    canal) : les protubérances sont des excès localisés au-dessus de cette
    ligne ; là où il n'y en a pas, rien n'est ajouté.

    `ha_sel` (0→1) pondère en plus ce détail par la ROUGEUR Hα : à 1, seules
    les structures rouges sont rehaussées, la couronne blanche est ignorée —
    c'est ce qui fait ressortir les protubérances là où la couronne est
    brillante et les noyait."""
    H, W = (base.shape[-2], base.shape[-1])
    yy, xx = np.indices((H, W)).astype(np.float32)
    r = np.hypot(xx - cx, yy - cy)
    w = smoothstep(R - 4.0, R + 2.0, r) * (1.0 - smoothstep(R + thickness,
                                                            R + thickness + 8.0, r))
    r_lo, r_hi = max(R - 4.0, 0.0), R + thickness + 8.0
    ring = (r >= r_lo) & (r <= r_hi)
    if not ring.any():
        return base.astype(np.float32), 0.0
    chans = [short_rad] if short_rad.ndim == 2 else \
            [short_rad[i] for i in range(short_rad.shape[0])]
    ri = r.astype(np.int32)
    ri_sel = ri[ring]
    bins = range(int(r_lo), int(r_hi) + 1)
    detail = []
    for c in chans:
        vals = c[ring]
        base_line = np.zeros(int(r_hi) + 2, np.float32)
        for b in bins:
            m = ri_sel == b
            if m.any():
                base_line[b] = float(np.median(vals[m]))   # attendu azimutal
        d = np.zeros((H, W), np.float32)
        d[ring] = np.maximum(vals - base_line[ri_sel], 0.0)  # excès = protubérance
        detail.append(d * w)
    det = detail[0] if len(detail) == 1 else np.stack(detail)
    # Sélection chromatique Hα : ne garde que ce qui est réellement rouge.
    if ha_sel > 1e-6 and short_rad.ndim == 3:
        hw = halpha_weight(short_rad, ring)
        mix = (1.0 - ha_sel) + ha_sel * hw          # ha_sel=0 → inchangé
        det = det * (mix if det.ndim == 2 else mix[None, :, :])
    # Purifier AVANT de saturer : le renfut ci-dessous amplifie la teinte
    # présente. Sur un détail déjà un peu jaune (V ≈ R), il écrase le bleu et
    # laisse V ≈ R → il FABRIQUE du jaune franc et baveux. En retirant d'abord
    # la fuite de vert, on ne sature plus qu'un rouge.
    if purify > 1e-6 and det.ndim == 3:
        det = purify_red(det, purify)
    # Renfort de saturation couleur : l'étirement final (GHS) écrase la
    # chrominance près du blanc ; pré-amplifier la couleur fait survivre le
    # rose des éruptions au stretch. Version MULTIPLICATIVE (puissance sur
    # les ratios canal/luminance) : saturation accrue en douceur, teinte
    # préservée, aucun canal annulé. chroma=1 → inchangé.
    if det.ndim == 3 and abs(chroma - 1.0) > 1e-6:
        dl = np.maximum(to_lum(det), 1e-12)[None, :, :]
        det = (dl * (np.maximum(det, 0.0) / dl) ** float(chroma)).astype(np.float32)
    det_lum = to_lum(det)
    pos = det_lum[det_lum > 0]
    ref_d = float(np.percentile(pos, 99.5)) if pos.size > 50 else 0.0
    ref_b = float(np.percentile(to_lum(base), 99.9))
    gain = (level_frac * ref_b / ref_d) if ref_d > 1e-12 else 0.0
    lay = det * gain
    if base.ndim == 3 and lay.ndim == 2:
        lay = np.stack([lay] * base.shape[0])
    if base.ndim == 2 and lay.ndim == 3:
        lay = to_lum(lay)
    return (base + lay).astype(np.float32), gain


def detect_stars(lum, cx, cy, r_min, k_sigma=6.0, max_stars=400):
    """Sources ponctuelles hors couronne : DoG, seuil kσ, maxima locaux,
    forme PSF (taille + rondeur). Retourne [(y, x)] triés par flux décroissant."""
    a = lum.astype(np.float32)
    H, W = a.shape
    yy, xx = np.indices((H, W)).astype(np.float32)
    r = np.hypot(xx - cx, yy - cy)
    zone = (r >= r_min)
    zone[:8, :] = False; zone[-8:, :] = False
    zone[:, :8] = False; zone[:, -8:] = False
    dog = gaussian_filter(a, 1.0) - gaussian_filter(a, 3.0)
    v = dog[zone]
    if v.size < 100:
        return []
    med = float(np.median(v))
    sig = 1.4826 * float(np.median(np.abs(v - med))) or 1e-12
    peaks = (dog == maximum_filter(dog, size=5)) & (dog > med + k_sigma * sig) & zone
    ys, xs = np.where(peaks)
    if len(ys) == 0:
        return []
    order = np.argsort(dog[ys, xs])[::-1][: max_stars * 3]
    stars = []
    for i in order:
        y, x = int(ys[i]), int(xs[i])
        if y < 6 or y >= H - 6 or x < 6 or x >= W - 6:
            continue
        p = a[y - 5:y + 6, x - 5:x + 6].astype(np.float64)
        bgv = float(np.median(np.concatenate([p[0, :], p[-1, :], p[:, 0], p[:, -1]])))
        q = np.maximum(p - bgv, 0.0)
        tot = q.sum()
        if tot <= 0:
            continue
        gy, gx = np.indices(q.shape)
        my = (q * gy).sum() / tot; mx = (q * gx).sum() / tot
        vy = (q * (gy - my) ** 2).sum() / tot
        vx = (q * (gx - mx) ** 2).sum() / tot
        sy, sx = np.sqrt(max(vy, 1e-9)), np.sqrt(max(vx, 1e-9))
        size = 0.5 * (sy + sx)
        round_ = abs(sy - sx) / max(sy, sx)
        if 0.4 <= size <= 3.0 and round_ <= 0.5:      # ponctuel et rond = étoile
            stars.append((y, x))
        if len(stars) >= max_stars:
            break
    return stars


def match_stars(list_a, list_b, tol=3.0):
    """Confirmation croisée : garde les détections présentes dans les deux
    listes à `tol` px près (élimine les pics de bruit)."""
    if not list_a or not list_b:
        return []
    b = np.asarray(list_b, np.float32)
    keep = []
    for (y, x) in list_a:
        d = np.hypot(b[:, 0] - y, b[:, 1] - x)
        if float(d.min()) <= tol:
            keep.append((y, x))
    return keep


def boost_stars(base, stack_rad, stars, level_frac):
    """Ajoute, à chaque étoile détectée, son motif PSF (fond local soustrait,
    fenêtré gaussien) pris dans l'empilement, amené à level_frac × P99,9(base).
    Seules les positions détectées sont touchées — le bruit n'est pas amplifié."""
    out = base.copy()
    lum = to_lum(stack_rad)
    H, W = lum.shape
    ref_b = float(np.percentile(to_lum(base), 99.9))
    gy, gx = np.indices((13, 13))
    win = np.exp(-(((gy - 6) ** 2 + (gx - 6) ** 2) / (2 * 2.2 ** 2))).astype(np.float32)
    # niveau : le pic médian des étoiles détectées est amené à level_frac×ref_b
    peaks = []
    for (y, x) in stars:
        if 6 <= y < H - 6 and 6 <= x < W - 6:
            p = lum[y - 6:y + 7, x - 6:x + 7]
            bgv = float(np.median(np.concatenate([p[0, :], p[-1, :], p[:, 0], p[:, -1]])))
            peaks.append(max(float(p[6, 6]) - bgv, 0.0))
    med_peak = float(np.median([p for p in peaks if p > 0] or [0.0]))
    if med_peak <= 1e-12:
        return out, 0.0
    gain = level_frac * ref_b / med_peak
    for (y, x) in stars:
        if not (6 <= y < H - 6 and 6 <= x < W - 6):
            continue
        p = lum[y - 6:y + 7, x - 6:x + 7]
        bgv = float(np.median(np.concatenate([p[0, :], p[-1, :], p[:, 0], p[:, -1]])))
        patch = np.maximum(p - bgv, 0.0) * win * gain
        if out.ndim == 2:
            out[y - 6:y + 7, x - 6:x + 7] += patch
        else:
            out[:, y - 6:y + 7, x - 6:x + 7] += patch[None, :, :]
    return out, gain


# ══════════════════════════════════════════════════════════════════════════════
#  Workers
# ══════════════════════════════════════════════════════════════════════════════
class BaseLoadWorker(QThread):
    done  = pyqtSignal(object, int, int)
    error = pyqtSignal(str)

    def __init__(self, siril, path=None):
        super().__init__()
        self.siril = siril
        self.path = path

    def run(self):
        tmp_path = None
        try:
            if self.path:
                src = self.path
            else:
                tmp = tempfile.mktemp(prefix="ecomp_")
                tmp_path = src = tmp + ".fit"
                with self.siril.image_lock():
                    self.siril.cmd(f'save "{tmp.replace(os.sep, "/")}"')
            raw, _ = load_fits(src)
            mx = float(np.max(raw))
            if mx > 1.0:
                raw = raw / (65535.0 if mx <= 65535.0 else mx)
            if raw.ndim == 3 and np.array_equal(raw[0], raw[1]) and \
               np.array_equal(raw[1], raw[2]):
                raw = raw[0]
            self.done.emit(raw, raw.shape[-1], raw.shape[-2])
        except Exception as e:
            self.error.emit(str(e))
        finally:
            if tmp_path and os.path.exists(tmp_path):
                try: os.remove(tmp_path)
                except Exception: pass


class ComposeWorker(QThread):
    progress = pyqtSignal(int, str)
    done     = pyqtSignal(object, object)      # résultat, étoiles [(y,x)]
    error    = pyqtSignal(str)

    def __init__(self, params, parent=None):
        super().__init__(parent)
        self.p = params

    def run(self):
        try:
            p = self.p
            base = p["base"].copy()
            cx, cy, R = p["cx"], p["cy"], p["R"]
            items = sorted(p["items"], key=lambda it: it[1])   # court → long
            stars_found = []

            stack_long = None
            if (p["moon_on"] or p["stars_on"]) and items:
                n = min(int(p["n_long"]), len(items))
                sel = items[-n:]
                self.progress.emit(10, f"Empilement des {n} poses longues "
                                       f"({', '.join(f'{e:g}s' for _, e in sel)})…")
                stack_long = median_stack_radiance(sel)
                if to_lum(stack_long).shape != (base.shape[-2], base.shape[-1]):
                    raise ValueError("poses et HDR de dimensions différentes")

            if p["moon_on"] and stack_long is not None:
                self.progress.emit(35, "Couche Lune — retrait du halo interne…")
                # Y a-t-il vraiment de la lumière cendrée enregistrée ? On compare
                # la structure DANS le disque au bruit : sans signal, amplifier ne
                # fait que remonter du grain — mieux vaut le dire franchement.
                lum_l = to_lum(stack_long)
                yy0, xx0 = np.indices(lum_l.shape).astype(np.float32)
                rr0 = np.hypot(xx0 - cx, yy0 - cy)
                inside = rr0 < 0.8 * R
                skip = False
                if inside.any():
                    hi = lum_l - gaussian_filter(lum_l, 2.0)
                    sig = 1.4826 * float(np.median(np.abs(hi - np.median(hi)))) or 1e-12
                    # La structure est mesurée sur l'image LISSÉE : il faut donc
                    # la comparer au bruit APRÈS ce même lissage, pas au bruit
                    # brut — sinon on rejette de la vraie lumière cendrée.
                    s_sm = 4.0
                    red = max(2.0 * s_sm * np.sqrt(np.pi), 1.0)   # réduction du bruit
                    struct = float(np.std(gaussian_filter(lum_l, s_sm)[inside]))
                    snr = struct / max(sig / red, 1e-12)
                    med_in = float(np.median(lum_l[inside]))
                    vide = med_in <= 2.0 * sig          # disque littéralement noir
                    if snr < 3.0 or vide:
                        raison = ("disque noir (aucun signal)" if vide
                                  else f"aucune structure au-dessus du bruit (SNR {snr:.2f})")
                        if p.get("moon_skip_dark", True):
                            skip = True
                            self.progress.emit(36,
                                f"Couche Lune IGNORÉE : {raison}. Aucune lumière cendrée "
                                f"n'a été enregistrée — l'amplifier n'ajouterait que du "
                                f"grain. Utilisez « Plancher du disque » pour un rendu.")
                        else:
                            self.progress.emit(36,
                                f"⚠ Disque lunaire : {raison}. La couche va amplifier "
                                f"du bruit — cochez « Ignorer si disque vide ».")
                    else:
                        self.progress.emit(36, f"Disque lunaire : structure détectée "
                                               f"(SNR {snr:.1f}) — clair de Terre exploitable.")
                if not skip:
                    layer = remove_inner_halo(stack_long, cx, cy, R) \
                        if p["moon_dehalo"] else stack_long
                    base, g = blend_moon(base, layer, cx, cy, R, p["moon_level"])
                    self.progress.emit(45, f"Couche Lune injectée (niveau "
                                           f"{100*p['moon_level']:.1f} %, gain ×{g:.3g}).")
            # Plancher uniforme : relève le disque à un gris sombre constant.
            # Purement cosmétique et assumé comme tel — cela n'invente aucun
            # détail, mais le disque cesse de lire comme un trou dans l'image.
            if p.get("moon_floor", 0.0) > 1e-6:
                yy0, xx0 = np.indices(base.shape[-2:]).astype(np.float32)
                rr0 = np.hypot(xx0 - cx, yy0 - cy)
                w = 1.0 - smoothstep(max(R - 10.0, 0.0), max(R - 2.0, 1.0), rr0)
                lvl = p["moon_floor"] * float(np.percentile(to_lum(base), 99.9))
                fl = (w * lvl).astype(np.float32)
                base = np.maximum(base, fl if base.ndim == 2 else fl[None, :, :])
                self.progress.emit(47, f"Plancher du disque : {lvl:.4g} "
                                       f"({100*p['moon_floor']:.2f} % de la couronne).")

            if p["prom_on"] and items:
                n = min(int(p["n_short"]), len(items))
                sel = items[:n]
                self.progress.emit(55, f"Couche protubérances — {n} pose(s) courte(s)…")
                short = median_stack_radiance(sel)
                base, g = blend_prominences(base, short, cx, cy, R,
                                            float(p["prom_thick"]), p["prom_level"],
                                            p.get("prom_chroma", 1.0),
                                            p.get("prom_ha", 0.0),
                                            p.get("purify", 0.0))
                self.progress.emit(65, f"Protubérances injectées (niveau "
                                       f"{100*p['prom_level']:.0f} %, gain ×{g:.3g}).")

            if p["stars_on"] and stack_long is not None:
                self.progress.emit(75, "Détection des étoiles (hors couronne)…")
                lum = to_lum(stack_long)
                stars = detect_stars(lum, cx, cy, p["star_rmin"], p["star_ksig"])
                if p["star_confirm"] and len(items) >= 2:
                    n = min(int(p["n_long"]), len(items))
                    sel = items[-n:]
                    half = max(1, len(sel) // 2)
                    sa = detect_stars(to_lum(median_stack_radiance(sel[:half])),
                                      cx, cy, p["star_rmin"], p["star_ksig"])
                    sb = detect_stars(to_lum(median_stack_radiance(sel[half:])),
                                      cx, cy, p["star_rmin"], p["star_ksig"])
                    confirmed = match_stars(sa, sb)
                    stars = match_stars(stars, confirmed, tol=3.0) if confirmed else []
                self.progress.emit(85, f"{len(stars)} étoile(s) retenue(s).")
                if stars:
                    base, g = boost_stars(base, stack_long, stars, p["star_level"])
                    self.progress.emit(90, f"Étoiles rehaussées (gain ×{g:.3g}).")
                stars_found = stars

            if (base.ndim == 3 and (p.get("neutral", 0.0) > 1e-6
                                     or p.get("neutral_wb", False))):
                self.progress.emit(94, "Neutralisation de la couronne…")
                base, g = neutralize_corona(base, cx, cy, R, p["neutral"],
                                            p.get("neutral_wb", True))
                self.progress.emit(95, f"Couronne neutralisée "
                                       f"(balance R/V/B ×{g[0]:.3f}/{g[1]:.3f}/{g[2]:.3f}).")

            if p.get("purify", 0.0) > 1e-6 and base.ndim == 3:
                before = base.copy()
                base = purify_red(base, p["purify"])
                touched = int((np.abs(before[1] - base[1]) > 1e-6).sum())
                self.progress.emit(96, f"Rouge pur : vert retiré sur {touched} px "
                                       f"chauds (jaune/orange → rouge).")

            # ── Garantie finale : le disque ne doit jamais ressortir noir ─────
            # Quoi qu'aient fait les couches en amont (la couche Lune remplace
            # le disque par l'empilement des poses longues, noir s'il n'y a pas
            # de clair de Terre), on le ramène ici au niveau du fond de ciel.
            # Idempotent : si Corona l'a déjà fait, rien ne change.
            if p.get("disc_grey", True) and R > 4:
                yy0, xx0 = np.indices(base.shape[-2:]).astype(np.float32)
                rr0 = np.hypot(xx0 - cx, yy0 - cy)
                lum_b = to_lum(base)
                # Fond de ciel estimé sur TOUTE l'image, indépendamment du rayon :
                # un anneau défini à partir de R sort du cadre si le rayon saisi
                # est trop grand, et la référence devient alors introuvable.
                # Le DÉCILE des pixels positifs donne le plancher de ciel :
                # un percentile plus haut (25 %) retombe encore dans la
                # couronne, très étendue, et le disque ressortirait trop clair.
                sky = lum_b > 0
                inside = rr0 < R * 0.9
                if sky.any() and inside.any():
                    lvl = float(np.percentile(lum_b[sky], 10.0))
                    cur = float(np.median(lum_b[inside]))
                    # Toujours appliqué : np.maximum ne peut qu'ÉCLAIRCIR ce qui
                    # est sous le niveau du fond. Un pixel déjà plus clair (y
                    # compris de la couronne, si le rayon saisi est trop grand)
                    # n'est pas touché — l'ancienne condition de déclenchement
                    # empêchait justement l'application dans ce cas.
                    w = (1.0 - smoothstep(R - 8.0, R - 1.0, rr0)).astype(np.float32)
                    ww = w if base.ndim == 2 else w[None, :, :]
                    base = np.maximum(base, lvl * ww).astype(np.float32)
                    new = float(np.median(to_lum(base)[inside]))
                    self.progress.emit(97, f"Disque : {cur:.5g} → {new:.5g} "
                                           f"(fond de ciel {lvl:.5g}).")
                else:
                    self.progress.emit(97, "⚠ Disque : pas de zone de fond de ciel "
                                           "exploitable — niveau inchangé.")
            elif p.get("disc_grey", True):
                self.progress.emit(97, f"⚠ Disque : rayon {R:.0f} px trop petit — "
                                       f"mise au niveau du fond non appliquée.")

            # sortie [0,1] : échelle préservée (division par le max seulement)
            mx = float(np.max(base))
            if mx > 1.0:
                base = base / mx
            self.progress.emit(100, "Composition terminée.")
            self.done.emit(np.maximum(base, 0.0).astype(np.float32), stars_found)
        except Exception as e:
            self.error.emit(str(e))


# ── Vue (pan/zoom, clic-glissé = centre, Maj+molette = rayon) ─────────────────────
class CenterView(QGraphicsView):
    center_at    = pyqtSignal(float, float)
    radius_delta = pyqtSignal(int)

    def __init__(self, scene, parent=None):
        super().__init__(scene, parent)
        self.setTransformationAnchor(QGraphicsView.ViewportAnchor.AnchorUnderMouse)
        self.setResizeAnchor(QGraphicsView.ViewportAnchor.AnchorUnderMouse)
        self.setDragMode(QGraphicsView.DragMode.NoDrag)
        self.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform)
        self._pan = None
        self._dragging = False

    def wheelEvent(self, event):
        up = event.angleDelta().y() > 0
        if event.modifiers() & Qt.KeyboardModifier.ShiftModifier:
            self.radius_delta.emit(1 if up else -1)
        else:
            f = 1.15 if up else 1 / 1.15
            self.scale(f, f)

    def mousePressEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            self._dragging = True
            pt = self.mapToScene(event.position().toPoint())
            self.center_at.emit(pt.x(), pt.y()); event.accept()
        elif event.button() == Qt.MouseButton.MiddleButton:
            self._pan = event.position().toPoint()
            self.setCursor(Qt.CursorShape.ClosedHandCursor); event.accept()
        else:
            super().mousePressEvent(event)

    def mouseMoveEvent(self, event):
        if self._dragging:
            pt = self.mapToScene(event.position().toPoint())
            self.center_at.emit(pt.x(), pt.y()); event.accept()
        elif self._pan is not None:
            d = event.position().toPoint() - self._pan
            self._pan = event.position().toPoint()
            self.horizontalScrollBar().setValue(self.horizontalScrollBar().value() - d.x())
            self.verticalScrollBar().setValue(self.verticalScrollBar().value() - d.y())
            event.accept()
        else:
            super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton and self._dragging:
            self._dragging = False; event.accept()
        elif event.button() == Qt.MouseButton.MiddleButton:
            self._pan = None; self.setCursor(Qt.CursorShape.ArrowCursor); event.accept()
        else:
            super().mouseReleaseEvent(event)


# ══════════════════════════════════════════════════════════════════════════════
#  Fenêtre principale
# ══════════════════════════════════════════════════════════════════════════════
class CompositeWindow(QMainWindow):

    def __init__(self, siril):
        super().__init__()
        self.siril = siril
        self.setWindowTitle(f"EclipseComposite — v{VERSION} (expérimental)")
        self.setMinimumSize(1100, 700)
        self.setStyleSheet(DARK_SS)

        self.base = None
        self.result = None
        self.cx = self.cy = 0.0
        self._H = self._W = 0
        self._items = []            # [(path, exptime)]
        self._pixitem = None
        self._overlay = []
        self._stars = []
        self._last_dir = os.getcwd()
        self._loader = None
        self._worker = None

        self._build_ui()
        self._load_base()

    # =========================================================================
    def _build_ui(self):
        central = QWidget(); self.setCentralWidget(central)
        h = QHBoxLayout(central); h.setContentsMargins(6, 6, 6, 6); h.setSpacing(8)

        left = QWidget(); left.setFixedWidth(330)
        lv = QVBoxLayout(left); lv.setContentsMargins(0, 0, 0, 0); lv.setSpacing(6)

        # Base HDR
        grp_b = QGroupBox("Image de base (HDR fusionné)")
        gb = QHBoxLayout(grp_b)
        b_cur = QPushButton("⟳ Image Siril"); b_cur.clicked.connect(self._load_base)
        b_file = QPushButton("Fichier…"); b_file.clicked.connect(self._open_base)
        gb.addWidget(b_cur); gb.addWidget(b_file)
        lv.addWidget(grp_b)

        # Poses alignées
        grp_p = QGroupBox("Poses alignées (celles de FusionHDR)")
        gp = QVBoxLayout(grp_p)
        self.table = QTableWidget(0, 2)
        self.table.setHorizontalHeaderLabels(["Fichier", "Expo (s)"])
        self.table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        self.table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setMaximumHeight(140)
        gp.addWidget(self.table)
        row = QHBoxLayout()
        row.addWidget(QLabel("Ajouter :"))
        for txt, fn in (("FITS", self._add_files), ("Dossier", self._add_folder)):
            b = QPushButton(txt); b.clicked.connect(fn); row.addWidget(b)
        row.addStretch()
        b_clr = QPushButton("Tout retirer"); b_clr.clicked.connect(self._clear_items)
        row.addWidget(b_clr)
        gp.addLayout(row)
        lv.addWidget(grp_p)

        # Géométrie
        grp_g = QGroupBox("Disque lunaire")
        gg = QGridLayout(grp_g)
        b_auto = QPushButton("◎ Détecter le centre"); b_auto.clicked.connect(self._auto_center)
        gg.addWidget(b_auto, 0, 0, 1, 2)
        self.lbl_center = QLabel("centre : —")
        self.lbl_center.setStyleSheet("font-size:9pt; color:#aaa;")
        gg.addWidget(self.lbl_center, 1, 0, 1, 2)
        gg.addWidget(QLabel("Rayon (px) :"), 2, 0)
        self.spn_r = QSpinBox(); self.spn_r.setRange(0, 20000)
        self.spn_r.valueChanged.connect(self._draw_overlay)
        gg.addWidget(self.spn_r, 2, 1)
        lv.addWidget(grp_g)

        # Couche Lune
        self.grp_moon = QGroupBox("☾ Lune — clair de Terre")
        self.grp_moon.setCheckable(True); self.grp_moon.setChecked(True)
        gm = QGridLayout(self.grp_moon)
        gm.addWidget(QLabel("Poses longues (N) :"), 0, 0)
        self.spn_nlong = QSpinBox(); self.spn_nlong.setRange(1, 50); self.spn_nlong.setValue(3)
        gm.addWidget(self.spn_nlong, 0, 1)
        self.chk_dehalo = QCheckBox("Retirer le halo interne")
        self.chk_dehalo.setChecked(True)
        gm.addWidget(self.chk_dehalo, 1, 0, 1, 2)
        gm.addWidget(QLabel("Niveau (% couronne) :"), 2, 0)
        self.spn_mlevel = QDoubleSpinBox(); self.spn_mlevel.setRange(0.1, 50.0)
        self.spn_mlevel.setSingleStep(0.5); self.spn_mlevel.setValue(3.0)
        self.spn_mlevel.setToolTip("La médiane du disque cendré est amenée à ce % de la "
                                   "couronne interne (P99,9 du HDR). Esthétique : ajustez.")
        gm.addWidget(self.spn_mlevel, 2, 1)
        gm.addWidget(QLabel("Plancher du disque (%) :"), 3, 0)
        self.spn_mfloor = QDoubleSpinBox(); self.spn_mfloor.setRange(0.0, 5.0)
        self.spn_mfloor.setSingleStep(0.05); self.spn_mfloor.setDecimals(2)
        self.spn_mfloor.setValue(0.0)
        self.spn_mfloor.setToolTip(
            "Relève le disque à un gris sombre uniforme, en % de la couronne.\n"
            "N'invente aucun détail — purement cosmétique — mais le disque cesse\n"
            "de lire comme un trou noir dans l'image. À utiliser quand vos poses\n"
            "longues ne contiennent PAS de lumière cendrée (le journal le dit).\n"
            "0,10 à 0,30 % suffit généralement.")
        gm.addWidget(self.spn_mfloor, 3, 1)
        self.chk_skip_dark = QCheckBox("Ignorer si disque vide (pas de clair de Terre)")
        self.chk_skip_dark.setChecked(True)
        self.chk_skip_dark.setToolTip(
            "Mesure la structure présente dans le disque et la compare au bruit.\n"
            "Si les poses longues y sont noires — cas courant : le clair de Terre\n"
            "demande plusieurs secondes de pose — la couche est simplement passée,\n"
            "au lieu d'amplifier du grain. Le journal indique la raison.\n"
            "Le « Plancher du disque » reste appliqué, lui.")
        gm.addWidget(self.chk_skip_dark, 4, 0, 1, 2)
        lv.addWidget(self.grp_moon)

        # Couche Protubérances
        self.grp_prom = QGroupBox("🔥 Protubérances — poses courtes")
        self.grp_prom.setCheckable(True); self.grp_prom.setChecked(True)
        gr = QGridLayout(self.grp_prom)
        gr.addWidget(QLabel("Poses courtes (N) :"), 0, 0)
        self.spn_nshort = QSpinBox(); self.spn_nshort.setRange(1, 20); self.spn_nshort.setValue(1)
        gr.addWidget(self.spn_nshort, 0, 1)
        gr.addWidget(QLabel("Épaisseur anneau (px) :"), 1, 0)
        self.spn_thick = QSpinBox(); self.spn_thick.setRange(2, 2000); self.spn_thick.setValue(4)
        self.spn_thick.setToolTip("Fin (2–8 px) = protubérances du limbe seulement ; "
                                  "large = risque de booster les racines des streamers.")
        gr.addWidget(self.spn_thick, 1, 1)
        gr.addWidget(QLabel("Niveau (% couronne) :"), 2, 0)
        self.spn_plevel = QDoubleSpinBox(); self.spn_plevel.setRange(1.0, 200.0)
        self.spn_plevel.setSingleStep(5.0); self.spn_plevel.setValue(80.0)
        gr.addWidget(self.spn_plevel, 2, 1)
        gr.addWidget(QLabel("Saturation couleur (×) :"), 3, 0)
        self.spn_chroma = QDoubleSpinBox(); self.spn_chroma.setRange(1.0, 5.0)
        self.spn_chroma.setSingleStep(0.5); self.spn_chroma.setValue(2.0)
        self.spn_chroma.setToolTip("Pré-amplifie la couleur des protubérances pour "
                                   "qu'elle survive à l'étirement GHS (qui écrase la "
                                   "chrominance près du blanc).")
        gr.addWidget(self.spn_chroma, 3, 1)
        gr.addWidget(QLabel("Sélectivité Hα (rouge) :"), 4, 0)
        self.spn_ha = QDoubleSpinBox(); self.spn_ha.setRange(0.0, 1.0)
        self.spn_ha.setSingleStep(0.1); self.spn_ha.setValue(0.7)
        self.spn_ha.setToolTip(
            "Ne rehausse que ce qui est réellement ROUGE (Hα). La couronne est\n"
            "blanche, les protubérances rouges : c'est le discriminant le plus\n"
            "sûr là où elles se noient dans une couronne brillante.\n"
            "0 = tout excès local ; 1 = uniquement le rouge.")
        gr.addWidget(self.spn_ha, 4, 1)
        lv.addWidget(self.grp_prom)

        # Couche Étoiles
        self.grp_star = QGroupBox("✶ Étoiles — détection + rehaussement")
        self.grp_star.setCheckable(True); self.grp_star.setChecked(False)
        gs = QGridLayout(self.grp_star)
        gs.addWidget(QLabel("Seuil (σ) :"), 0, 0)
        self.spn_ksig = QDoubleSpinBox(); self.spn_ksig.setRange(3.0, 20.0)
        self.spn_ksig.setSingleStep(0.5); self.spn_ksig.setValue(6.0)
        gs.addWidget(self.spn_ksig, 0, 1)
        gs.addWidget(QLabel("Rayon min (× R lune) :"), 1, 0)
        self.spn_rmin = QDoubleSpinBox(); self.spn_rmin.setRange(1.2, 20.0)
        self.spn_rmin.setSingleStep(0.1); self.spn_rmin.setValue(2.5)
        self.spn_rmin.setToolTip("En deçà, la couronne domine : pas de détection.")
        gs.addWidget(self.spn_rmin, 1, 1)
        self.chk_confirm = QCheckBox("Confirmation croisée (2 demi-empilements)")
        self.chk_confirm.setChecked(True)
        gs.addWidget(self.chk_confirm, 2, 0, 1, 2)
        gs.addWidget(QLabel("Niveau (% couronne) :"), 3, 0)
        self.spn_slevel = QDoubleSpinBox(); self.spn_slevel.setRange(1.0, 100.0)
        self.spn_slevel.setSingleStep(5.0); self.spn_slevel.setValue(25.0)
        gs.addWidget(self.spn_slevel, 3, 1)
        lv.addWidget(self.grp_star)

        grp_neu = QGroupBox("Couronne — couleur")
        gn = QGridLayout(grp_neu)
        gn.addWidget(QLabel("Neutraliser (Hα préservé) :"), 0, 0)
        self.spn_neutral = QDoubleSpinBox(); self.spn_neutral.setRange(0.0, 1.0)
        self.spn_neutral.setSingleStep(0.1); self.spn_neutral.setValue(0.8)
        self.spn_neutral.setToolTip(
            "La couronne K est physiquement BLANCHE : toute teinte est un\n"
            "artefact. La désaturer maximise le RAPPORT de saturation avec les\n"
            "protubérances rouges — la saturation absolue se remonte ensuite\n"
            "d'un curseur dans n'importe quel logiciel photo.\n"
            "0 = inchangé · 1 = couronne parfaitement neutre.")
        gn.addWidget(self.spn_neutral, 0, 1)
        self.chk_wb = QCheckBox("Retirer la dominante (balance sur la couronne)")
        self.chk_wb.setChecked(True)
        self.chk_wb.setToolTip("Égalise les canaux sur la couronne seule, en excluant "
                               "les pixels rouges pour ne pas blanchir les protubérances.")
        gn.addWidget(self.chk_wb, 1, 0, 1, 2)
        gn.addWidget(QLabel("Rouge pur (retire le vert) :"), 2, 0)
        self.spn_purify = QDoubleSpinBox(); self.spn_purify.setRange(0.0, 1.0)
        self.spn_purify.setSingleStep(0.1); self.spn_purify.setValue(0.0)
        self.spn_purify.setToolTip(
            "Sur les pixels où le rouge domine le bleu, ramène le vert vers le\n"
            "bleu : le jaune/orange redevient un rouge franc. Le blanc et le\n"
            "bleu ne bougent pas.\n"
            "Plus agressif que le SCNR de Siril, qui se contente de plafonner\n"
            "le vert au niveau neutre. 0 = inchangé.")
        gn.addWidget(self.spn_purify, 2, 1)
        self.chk_disc = QCheckBox("Disque au niveau du fond (jamais noir)")
        self.chk_disc.setChecked(True)
        self.chk_disc.setToolTip(
            "Garantie finale : si une couche a laissé le disque plus sombre que\n"
            "le fond de ciel (la couche Lune le remplace par l'empilement des\n"
            "poses longues, noir s'il n'y a pas de clair de Terre), il est\n"
            "ramené au niveau du fond. Sans effet si Corona l'a déjà fait.")
        gn.addWidget(self.chk_disc, 3, 0, 1, 2)
        lv.addWidget(grp_neu)

        hint = QLabel("Clic-glissé = centre · Maj+molette = rayon · molette = zoom. "
                      "Composite esthétique : niveaux à ajuster à l'œil, puis GHS dans Siril.")
        hint.setWordWrap(True); hint.setStyleSheet("font-size:8pt; color:#888;")
        lv.addWidget(hint)

        lv.addStretch()
        self.btn_go = QPushButton("▶  Composer"); self.btn_go.setObjectName("BtnGo")
        self.btn_go.clicked.connect(self._compose)
        lv.addWidget(self.btn_go)
        self.btn_save = QPushButton("Enregistrer + charger dans Siril")
        self.btn_save.setEnabled(False); self.btn_save.clicked.connect(self._save_load)
        lv.addWidget(self.btn_save)
        h.addWidget(left)

        right = QWidget()
        rv = QVBoxLayout(right); rv.setContentsMargins(0, 0, 0, 0); rv.setSpacing(4)
        self._scene = QGraphicsScene()
        self._view = CenterView(self._scene)
        self._view.center_at.connect(self._on_center_at)
        self._view.radius_delta.connect(self._on_radius_delta)
        rv.addWidget(self._view, stretch=1)
        self.prog = QProgressBar(); rv.addWidget(self.prog)
        self.log = QTextEdit(); self.log.setReadOnly(True); self.log.setFixedHeight(100)
        rv.addWidget(self.log)
        h.addWidget(right, stretch=1)

        self.statusBar().showMessage("Chargement…")

    # =========================================================================
    #  Base & poses
    # =========================================================================
    def _load_base(self, path=None):
        self.statusBar().showMessage("Chargement de l'image de base…")
        self._loader = BaseLoadWorker(self.siril, path if isinstance(path, str) else None)
        self._loader.done.connect(self._on_base)
        self._loader.error.connect(lambda e: self._log(f"✗ {e}"))
        self._loader.start()

    def _open_base(self):
        f, _ = QFileDialog.getOpenFileName(self, "Ouvrir le HDR", os.getcwd(),
                                           "FITS (*.fit *.fits *.fts *.FIT *.FITS *.FTS)")
        if f:
            self._load_base(f)

    def _on_base(self, data, W, H):
        self.base = data
        self.result = None
        self._H, self._W = H, W
        self.btn_save.setEnabled(False)
        self._show(autostretch_display(data))
        self._log(f"Base {W}×{H} ({'couleur' if data.ndim == 3 else 'mono'}) chargée.")
        self._auto_center()

    def _add_paths(self, paths):
        added = 0
        for p in paths:
            try:
                hdr = fits.getheader(p)
                e = get_exptime(hdr)
                self._items.append((p, e))
                added += 1
            except Exception as ex:
                self._log(f"✗ {Path(p).name} : {ex}")
        if added:
            self._items.sort(key=lambda it: it[1])
            self._refresh_table()
            miss = sum(1 for _, e in self._items if e <= 0)
            if miss:
                self._log(f"⚠ {miss} pose(s) sans EXPTIME — elles seront traitées "
                          "comme les plus courtes.")
            self._log(f"+ {added} pose(s), triées par exposition.")

    def _add_files(self):
        files, _ = QFileDialog.getOpenFileNames(self, "Ajouter des poses",
            self._last_dir, "FITS (*.fit *.fits *.fts *.FIT *.FITS *.FTS)")
        if files:
            self._last_dir = str(Path(files[0]).parent)
            self._add_paths(files)

    def _add_folder(self):
        d = QFileDialog.getExistingDirectory(self, "Dossier des poses", self._last_dir)
        if not d:
            return
        self._last_dir = d
        files = sorted(str(f) for f in Path(d).iterdir()
                       if f.is_file() and f.suffix.lower() in FITS_EXT)
        if files:
            self._add_paths(files)
        else:
            self._log("Aucun FITS dans ce dossier.")

    def _clear_items(self):
        self._items.clear()
        self._refresh_table()

    def _refresh_table(self):
        self.table.setRowCount(len(self._items))
        for i, (p, e) in enumerate(self._items):
            nm = QTableWidgetItem(Path(p).name); nm.setToolTip(p)
            ex = QTableWidgetItem(f"{e:g}" if e > 0 else "—")
            for c, itm in enumerate((nm, ex)):
                itm.setFlags(itm.flags() & ~Qt.ItemFlag.ItemIsEditable)
                self.table.setItem(i, c, itm)

    # =========================================================================
    #  Géométrie
    # =========================================================================
    def _auto_center(self):
        if self.base is None:
            return
        d, src = None, ""

        # 1) Cercle écrit par SirilJ Align dans l'en-tête des poses : c'est la
        #    source la plus sûre, mesurée sur la donnée brute au moment du
        #    recalage. (Absent des fichiers alignés avant cette version.)
        for path, _e in self._items:
            try:
                h = fits.getheader(path)
                if all(k in h for k in ("MOONX", "MOONY", "MOONR")):
                    d = (float(h["MOONX"]), float(h["MOONY"]), float(h["MOONR"]))
                    src = f"en-tête de {Path(path).name}"
                    break
            except Exception:
                continue

        # 2) Détection sur les POSES, jamais sur l'image de base : celle-ci est
        #    déjà traitée (Corona), et un « rayon max couronne » y laisse un
        #    cercle parfait — bien plus net que le limbe — sur lequel
        #    l'ajustement se fixe, d'où des rayons absurdes.
        #    Consensus par MÉDIANE sur plusieurs poses : une détection isolée
        #    aberrante est ainsi écartée d'office.
        if d is None and self._items:
            self._log("Recherche du disque sur les poses (source la plus fiable)…")
            found = []
            for path, _e in sorted(self._items, key=lambda it: -it[1])[:6]:
                try:
                    data, _ = load_fits(path)
                    if data.shape[-2:] != (self._H, self._W):
                        continue
                    c = detect_disk(to_lum(data))
                    if c is not None:
                        found.append(c)
                except Exception:
                    continue
            if found:
                rs = np.array([c[2] for c in found])
                r_med = float(np.median(rs))
                keep = [c for c in found if abs(c[2] - r_med) < 0.15 * r_med]
                if not keep:
                    keep = found
                d = (float(np.median([c[0] for c in keep])),
                     float(np.median([c[1] for c in keep])),
                     float(np.median([c[2] for c in keep])))
                src = f"médiane de {len(keep)}/{len(found)} pose(s)"
                if len(found) > len(keep):
                    self._log(f"   {len(found)-len(keep)} détection(s) aberrante(s) écartée(s).")

        # 3) En dernier recours seulement : l'image de base.
        if d is None:
            c = detect_disk(to_lum(self.base))
            if c is not None:
                d, src = c, "image de base (à vérifier)"

        if d is not None:
            self.cx, self.cy, r = d
            self.spn_r.blockSignals(True)
            self.spn_r.setValue(int(round(r)))
            self.spn_r.blockSignals(False)
            self._log(f"Disque détecté ({src}) : centre ({self.cx:.0f}, {self.cy:.0f}), "
                      f"rayon {r:.0f} px.")
        else:
            self.cx, self.cy = self._W / 2.0, self._H / 2.0
            self._log("⚠ Disque non détecté, même sur les poses — centre au milieu. "
                      "Réglez au clic-glissé, et le rayon à Maj+molette.")
        self.lbl_center.setText(f"centre : ({self.cx:.0f}, {self.cy:.0f})")
        self._draw_overlay()

    def _on_center_at(self, sx, sy):
        if self.base is None:
            return
        self.cx = float(sx); self.cy = float(self._H - 1 - sy)
        self.lbl_center.setText(f"centre : ({self.cx:.0f}, {self.cy:.0f})")
        self._draw_overlay()

    def _on_radius_delta(self, direction):
        step = max(1, int(self.spn_r.value() * 0.03))
        self.spn_r.setValue(self.spn_r.value() + direction * step)

    def _draw_overlay(self):
        for it in self._overlay:
            self._scene.removeItem(it)
        self._overlay.clear()
        if self.base is None:
            return
        sx = self.cx; sy = self._H - 1 - self.cy
        pen = QPen(QColor("#ff3333")); pen.setCosmetic(True); pen.setWidth(2)
        L = 12
        self._overlay.append(self._scene.addLine(sx - L, sy, sx + L, sy, pen))
        self._overlay.append(self._scene.addLine(sx, sy - L, sx, sy + L, pen))
        R = self.spn_r.value()
        if R > 0:
            penc = QPen(QColor("#33aaff")); penc.setCosmetic(True); penc.setWidth(2)
            self._overlay.append(self._scene.addEllipse(sx - R, sy - R, 2*R, 2*R, penc))
        # étoiles retenues (après composition)
        pens = QPen(QColor("#ffe066")); pens.setCosmetic(True); pens.setWidth(1)
        for (y, x) in self._stars[:400]:
            syy = self._H - 1 - y
            self._overlay.append(self._scene.addEllipse(x - 6, syy - 6, 12, 12, pens))

    def _show(self, disp):
        pm = to_pixmap(disp)
        if self._pixitem is None:
            self._pixitem = self._scene.addPixmap(pm)
        else:
            self._pixitem.setPixmap(pm)
        self._scene.setSceneRect(0, 0, pm.width(), pm.height())
        QTimer.singleShot(50, lambda: self._view.fitInView(
            self._scene.sceneRect(), Qt.AspectRatioMode.KeepAspectRatio))
        self._draw_overlay()

    # =========================================================================
    #  Composition
    # =========================================================================
    def _compose(self):
        if self.base is None:
            self._log("Aucune image de base."); return
        R = float(self.spn_r.value())
        if R <= 0:
            self._log("⚠ Rayon lunaire nul — cliquez « Détecter le centre »."); return
        needs = (self.grp_moon.isChecked() or self.grp_prom.isChecked()
                 or self.grp_star.isChecked())
        if not needs:
            self._log("Aucune couche activée."); return
        if not self._items:
            self._log("⚠ Ajoutez les poses alignées (FITS ou Dossier)."); return

        params = {
            "base": self.base, "items": list(self._items),
            "cx": self.cx, "cy": self.cy, "R": R,
            "moon_on":  self.grp_moon.isChecked(),
            "n_long":   self.spn_nlong.value(),
            "moon_dehalo": self.chk_dehalo.isChecked(),
            "moon_level":  self.spn_mlevel.value() / 100.0,
            "moon_floor":  self.spn_mfloor.value() / 100.0,
            "moon_skip_dark": self.chk_skip_dark.isChecked(),
            "prom_on":  self.grp_prom.isChecked(),
            "n_short":  self.spn_nshort.value(),
            "prom_thick": float(self.spn_thick.value()),
            "prom_level": self.spn_plevel.value() / 100.0,
            "prom_chroma": self.spn_chroma.value(),
            "prom_ha":     self.spn_ha.value(),
            "neutral":     self.spn_neutral.value(),
            "neutral_wb":  self.chk_wb.isChecked(),
            "purify":      self.spn_purify.value(),
            "disc_grey":   self.chk_disc.isChecked(),
            "stars_on": self.grp_star.isChecked(),
            "star_ksig": self.spn_ksig.value(),
            "star_rmin": self.spn_rmin.value() * R,
            "star_confirm": self.chk_confirm.isChecked(),
            "star_level": self.spn_slevel.value() / 100.0,
        }
        self.btn_go.setEnabled(False)
        self.btn_save.setEnabled(False)
        self.prog.setValue(0)
        self._log("━━ Composition…")
        self._worker = ComposeWorker(params)
        self._worker.progress.connect(lambda v, m: (self.prog.setValue(v), self._log(m)))
        self._worker.done.connect(self._on_done)
        self._worker.error.connect(lambda e: (self._log(f"✗ {e}"),
                                              self.btn_go.setEnabled(True)))
        self._worker.start()

    def _on_done(self, result, stars):
        self.result = result
        self._stars = list(stars)
        self.btn_go.setEnabled(True)
        self.btn_save.setEnabled(True)
        self._show(autostretch_display(result))
        self._log("✓ Composition prête (aperçu autostretch). Ajustez les niveaux "
                  "ou enregistrez, puis GHS dans Siril.")

    def _save_load(self):
        if self.result is None:
            return
        default = str(Path(os.getcwd()) / "hdr_composite.fit")   # dossier Siril
        f, _ = QFileDialog.getSaveFileName(self, "Enregistrer le composite",
                                           default, "FITS (*.fit *.fits)")
        if not f:
            return
        try:
            save_fits(self.result, f)
            self._log(f"✓ Enregistré : {f}")
            try:
                try:
                    self.siril.cmd("set32bits")
                except Exception:
                    pass
                self.siril.cmd(f'load "{Path(f).with_suffix("")}"')
                self.siril.log(f"EclipseComposite — chargé : {Path(f).name}")
                self._log("✓ Chargé dans Siril (32 bits).")
            except Exception:
                self._log("✓ Enregistré (chargez-le manuellement).")
        except Exception as e:
            self._log(f"✗ {e}")

    # =========================================================================
    def _log(self, msg):
        self.log.append(msg)
        self.log.ensureCursorVisible()
        self.statusBar().showMessage(msg)

    def closeEvent(self, event):
        for w in (self._loader, self._worker):
            if w and w.isRunning():
                w.requestInterruption()
                w.wait(4000)
        super().closeEvent(event)


def main():
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    siril = s.SirilInterface()
    try:
        siril.connect()
    except Exception as e:
        print(f"Avertissement : connexion Siril impossible ({e}).")
    win = CompositeWindow(siril)
    win.show()
    app.exec()


if __name__ == "__main__":
    main()
