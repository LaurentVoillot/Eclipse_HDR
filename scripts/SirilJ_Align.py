##############################################
# SirilJ — Align (Soleil & Éclipse)
# Alignement unifié : surface solaire OU disque lunaire
# Version 1.0.0
##############################################
#
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Fusionne SirilJ_SolarAlign_v1.1 (multi-points, surface) et EclipseAlign
# (recalage sur la Lune, inter-poses) en un seul outil à deux modes. Les
# deux scripts d'origine restent inchangés.
#
# Crédits
# -------
#   • Soleil  : corrélation de phase multi-points (OpenCV), à la
#     AutoStakkert! (lucky imaging solaire)
#   • Éclipse : détection du limbe lunaire par ajustement de cercle
#     RANSAC + Kåsa, recalage par translation
#   • Interface : PyQt6 / conventions VeraLux

"""
SirilJ Align — un seul programme, deux cibles.

  ☼ Soleil (surface, multi-points)
    Aligne la séquence courante de Siril (même exposition) sur la structure
    de surface (granulation, protubérances). Translation globale par
    corrélation de phase + champ de déformation local multi-points qui
    corrige la turbulence. Sortie : séquence alignée → empilement.

  ☾ Éclipse (Lune, inter-poses)
    Recale des poses de luminosités TRÈS différentes (du diamant à la
    couronne) en les alignant sur le seul repère commun : le disque
    lunaire. Détection adaptative (centroïde du disque enclos ou limbe
    RANSAC) + translation. Sortie : poses recalées → FusionHDR.

Compatibilité
-------------
• Siril 1.3+   • Python 3.10+ (via sirilpy)
• Dépendances : numpy, opencv-python, astropy, scipy, PyQt6
"""

import sys
import os
import re
import html

try:
    import sirilpy as s
except ImportError:
    print("Erreur : module sirilpy introuvable.")
    sys.exit(1)

s.ensure_installed("numpy", "astropy", "scipy", "opencv-python", "PyQt6")

import numpy as np
import cv2
from pathlib import Path
from typing import Optional

from astropy.io import fits
from scipy.ndimage import gaussian_filter, shift as ndshift, label, binary_closing

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
    QLabel, QLineEdit, QPushButton, QComboBox, QCheckBox, QSpinBox,
    QProgressBar, QTextEdit, QFileDialog, QGroupBox, QStackedWidget,
    QGraphicsView, QGraphicsScene, QGraphicsPixmapItem,
    QTableWidget, QTableWidgetItem, QHeaderView, QAbstractItemView
)
from PyQt6.QtCore import Qt, QThread, pyqtSignal, QSettings, QTimer
from PyQt6.QtGui import QImage, QPixmap, QPainter, QPen, QColor

VERSION = "1.0.0"

DARK_SS = """
QWidget            { background:#2b2b2b; color:#d4d4d4; font-size:11px; }
QMainWindow        { background:#2b2b2b; }
QPushButton        { background:#3c3c3c; color:#d4d4d4; border:1px solid #555;
                     border-radius:4px; padding:4px 10px; }
QPushButton:hover  { background:#4a4a4a; }
QPushButton:pressed{ background:#555; }
QPushButton:disabled { color:#666; border-color:#444; }
QPushButton#BtnGo  { background:#1e6b1e; color:#fff; font-weight:bold;
                     border:1px solid #2a9b2a; padding:6px 12px; }
QPushButton#BtnGo:hover    { background:#248c24; }
QPushButton#BtnGo:disabled { background:#333; color:#666; border-color:#444; }
QPushButton#BtnStart { background:#1e6b1e; color:#fff; font-weight:bold;
                       font-size:12px; border:1px solid #2a9b2a; padding:6px 12px; }
QPushButton#BtnStart:hover    { background:#248c24; }
QPushButton#BtnStart:disabled { background:#333; color:#666; border-color:#444; }
QPushButton#BtnPrep  { background:#1e3d6b; color:#fff; border:1px solid #2a60a0;
                       padding:4px 10px; }
QPushButton#BtnPrep:hover    { background:#255080; }
QPushButton#BtnPrep:disabled { background:#333; color:#666; border-color:#444; }
QLineEdit          { background:#1e1e1e; border:1px solid #555; border-radius:3px;
                     padding:3px 6px; color:#d4d4d4; }
QComboBox          { background:#3c3c3c; border:1px solid #555; border-radius:3px;
                     padding:3px 6px; color:#d4d4d4; }
QComboBox QAbstractItemView { background:#2b2b2b; color:#d4d4d4;
                               selection-background-color:#555; }
QSpinBox           { background:#1e1e1e; border:1px solid #555; border-radius:3px;
                     padding:2px 4px; color:#d4d4d4; }
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
QScrollBar:horizontal { background:#2b2b2b; height:10px; }
QScrollBar::handle:horizontal { background:#555; border-radius:5px; }
QLabel#SeqInfoOK   { color:#7ec77e; font-weight:bold; }
QLabel#SeqInfoKO   { color:#d77; font-weight:bold; }
"""

FITS_EXT = {".fit", ".fits", ".fts"}


# ── FITS / display helpers ────────────────────────────────────────────────────
def natural_key(name: str):
    return [int(c) if c.isdigit() else c.lower() for c in re.split(r"(\d+)", name)]

def load_fits(path):
    """FITS → (data float32 aux valeurs d'origine, header). Couleur → (C,H,W)."""
    with fits.open(path) as hdul:
        hdr = hdul[0].header
        data = hdul[0].data.astype(np.float32)
    if data.ndim == 3 and data.shape[2] in (3, 4) and data.shape[0] not in (3, 4):
        data = data.transpose(2, 0, 1)
    if data.ndim == 3 and data.shape[0] == 4:
        data = data[:3]
    return data, hdr

def to_lum(d):
    return d if d.ndim == 2 else d.mean(axis=0)

def reg_lum(data):
    """Luminance normalisée [0,1] pour l'enregistrement / la corrélation."""
    lum = to_lum(data).astype(np.float32)
    mn, mx = float(lum.min()), float(lum.max())
    return (lum - mn) / (mx - mn) if mx > mn else lum

def estimate_sharpness(lum):
    return float(np.var(np.diff(lum, axis=0)) + np.var(np.diff(lum, axis=1)))

def save_fits(data, hdr, path):
    """FITS float32 (BITPIX=-32), header préservé (dont EXPTIME pour FusionHDR)."""
    out = np.asarray(data, dtype=np.float32)
    h = hdr.copy()
    for kw in ("BZERO", "BSCALE", "BLANK", "DATAMAX", "DATAMIN"):
        h.remove(kw, ignore_missing=True)
    fits.PrimaryHDU(out, header=h).writeto(path, overwrite=True)

def _mtf(x, m, lo, hi):
    dist = hi - lo
    if dist < 1e-9:
        return np.where(x > lo, 1.0, 0.0).astype(np.float32)
    xp = np.clip((x - lo) / dist, 0.0, 1.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        y = ((m - 1.0) * xp) / ((2.0 * m - 1.0) * xp - m)
    return np.nan_to_num(y, nan=0.0, posinf=1.0, neginf=0.0).astype(np.float32)

def autostretch_lum(lum):
    x = lum.astype(np.float32)
    mn, mx = float(x.min()), float(x.max())
    x = (x - mn) / (mx - mn) if mx > mn else np.zeros_like(x)
    st = max(1, x.size // 500_000)
    samp = x.flatten()[::st]
    med = float(np.median(samp))
    mad = float(np.median(np.abs(samp - med))) * 1.4826 or 0.001
    c0 = max(0.0, med - 2.8 * mad)
    mt = float(_mtf(np.float32(med - c0), 0.25, 0.0, 1.0))
    return np.clip(_mtf(x, mt, c0, 1.0), 0, 1)

def lum_to_pixmap(lum01):
    # FITS : ligne 0 = bas → flip vertical pour afficher comme Siril.
    g = np.ascontiguousarray(np.flipud((np.clip(lum01, 0, 1) * 255).astype(np.uint8)))
    h, w = g.shape
    return QPixmap.fromImage(QImage(g.tobytes(), w, h, w, QImage.Format.Format_Grayscale8))


# ══════════════════════════════════════════════════════════════════════════════
#  CŒUR SOLEIL — alignement multi-points (corrélation de phase)
# ══════════════════════════════════════════════════════════════════════════════
MIN_RESP = 0.03   # réponse minimale de la corrélation de phase (pic normalisé)


def noise_sigma(lum):
    """Bruit robuste : 1,4826 × MAD(image − flou gaussien σ=1)."""
    hi = lum - cv2.GaussianBlur(lum, (0, 0), 1.0)
    return 1.4826 * float(np.median(np.abs(hi - np.median(hi))))


def make_aps(ref, box, step):
    """Grille de points d'ancrage sur la référence (luminance [0,1]).

    Un point est retenu si son patch a une structure au-dessus du bruit
    (σ_patch > 2,5 × σ_bruit) → exclut le ciel vide, garde granulation,
    protubérances et limbe.
    """
    H, W = ref.shape
    sig = max(noise_sigma(ref), 1e-6)
    half = box // 2
    aps = []
    for cy in range(half, H - half + 1, step):
        for cx in range(half, W - half + 1, step):
            patch = ref[cy - half:cy + half, cx - half:cx + half]
            if float(patch.std()) > 2.5 * sig:
                aps.append((cy, cx))
    return aps


def _hann(shape_wh):
    return cv2.createHanningWindow(shape_wh, cv2.CV_32F)


def global_shift(ref, mov, hann):
    """Translation globale (dx, dy) de mov par rapport à ref (sub-pixel).

    ATTENTION : ne JAMAIS passer la fenêtre à cv2.phaseCorrelate — OpenCV
    multiplierait alors ses entrées EN PLACE (corruption silencieuse).
    On pré-multiplie hors-place : résultat identique, entrées intactes.
    """
    (dx, dy), resp = cv2.phaseCorrelate(ref * hann, mov * hann)
    return float(dx), float(dy), float(resp)


def ap_shifts(ref, mov, aps, box, gdx, gdy, max_local):
    """Décalage total (global + résiduel local) à chaque point d'ancrage.
    Retourne (dxs, dys, valid) alignés sur `aps`."""
    H, W = ref.shape
    half = box // 2
    ox, oy = int(round(gdx)), int(round(gdy))
    hann_p = _hann((box, box))
    n = len(aps)
    dxs = np.zeros(n, np.float32)
    dys = np.zeros(n, np.float32)
    valid = np.zeros(n, bool)
    for i, (cy, cx) in enumerate(aps):
        rp = ref[cy - half:cy + half, cx - half:cx + half]
        # patch homologue dans mov, décalé du global (borné dans l'image)
        my0 = min(max(cy + oy - half, 0), H - box)
        mx0 = min(max(cx + ox - half, 0), W - box)
        mp_ = mov[my0:my0 + box, mx0:mx0 + box]
        # pré-multiplication hors-place (jamais de fenêtre passée à OpenCV,
        # qui modifierait rp/mp_ en place — vues sur ref/mov !)
        (ldx, ldy), resp = cv2.phaseCorrelate(rp * hann_p, mp_ * hann_p)
        tdx = (mx0 + half - cx) + ldx
        tdy = (my0 + half - cy) + ldy
        if resp >= MIN_RESP and abs(tdx - gdx) <= max_local and abs(tdy - gdy) <= max_local:
            dxs[i], dys[i], valid[i] = tdx, tdy, True
    return dxs, dys, valid


def build_maps(shape, aps, dxs, dys, valid, gdx, gdy, grid_sigma=0.75):
    """Champ de déformation dense (map_x, map_y) pour cv2.remap.

    Décalages épars posés sur leur grille, invalides comblés par moyenne
    pondérée des voisins (flou masqué), grille lissée puis interpolée
    bilinéairement en pleine résolution (extrapolation constante aux bords).
    """
    H, W = shape
    ys = sorted({cy for cy, _ in aps})
    xs = sorted({cx for _, cx in aps})
    iy = {v: i for i, v in enumerate(ys)}
    ix = {v: i for i, v in enumerate(xs)}
    gh, gw = len(ys), len(xs)
    gdx_grid = np.zeros((gh, gw), np.float32)
    gdy_grid = np.zeros((gh, gw), np.float32)
    mask = np.zeros((gh, gw), np.float32)
    for (cy, cx), dx, dy, ok in zip(aps, dxs, dys, valid):
        if ok:
            gdx_grid[iy[cy], ix[cx]] = dx
            gdy_grid[iy[cy], ix[cx]] = dy
            mask[iy[cy], ix[cx]] = 1.0
    if mask.sum() == 0:
        gdx_grid[:] = gdx
        gdy_grid[:] = gdy
    else:
        sig = max(grid_sigma, 0.6)
        wm = cv2.GaussianBlur(mask, (0, 0), sig, borderType=cv2.BORDER_REPLICATE)
        wm = np.maximum(wm, 1e-6)
        gdx_grid = cv2.GaussianBlur(gdx_grid * mask, (0, 0), sig,
                                    borderType=cv2.BORDER_REPLICATE) / wm
        gdy_grid = cv2.GaussianBlur(gdy_grid * mask, (0, 0), sig,
                                    borderType=cv2.BORDER_REPLICATE) / wm

    ys_a = np.asarray(ys, np.float32)
    xs_a = np.asarray(xs, np.float32)
    yy = np.arange(H, dtype=np.float32)
    xx = np.arange(W, dtype=np.float32)
    fy = np.interp(yy, ys_a, np.arange(gh, dtype=np.float32))
    fx = np.interp(xx, xs_a, np.arange(gw, dtype=np.float32))
    y0 = np.clip(np.floor(fy).astype(np.int32), 0, gh - 1)
    x0 = np.clip(np.floor(fx).astype(np.int32), 0, gw - 1)
    y1 = np.minimum(y0 + 1, gh - 1)
    x1 = np.minimum(x0 + 1, gw - 1)
    wy = (fy - y0).astype(np.float32)[:, None]
    wx = (fx - x0).astype(np.float32)[None, :]

    def bilin(g):
        a = g[np.ix_(y0, x0)] * (1 - wy) * (1 - wx)
        b = g[np.ix_(y0, x1)] * (1 - wy) * wx
        c = g[np.ix_(y1, x0)] * wy * (1 - wx)
        d = g[np.ix_(y1, x1)] * wy * wx
        return (a + b + c + d).astype(np.float32)

    DX = bilin(gdx_grid)
    DY = bilin(gdy_grid)
    map_x = np.arange(W, dtype=np.float32)[None, :] + DX
    map_y = np.arange(H, dtype=np.float32)[:, None] + DY
    return map_x, map_y


def warp_channels(data, map_x, map_y):
    fn = lambda ch: cv2.remap(ch.astype(np.float32), map_x, map_y,
                              cv2.INTER_CUBIC, borderMode=cv2.BORDER_REPLICATE)
    if data.ndim == 2:
        return fn(data)
    return np.stack([fn(data[c]) for c in range(data.shape[0])])


def mp_maps(ref, mov, aps, box, max_local, hann, n_pass=2, grid_sigma=0.75):
    """Cartes de déformation complètes (global + local), raffinées en n_pass.

    Chaque passe mesure le résiduel sur la luminance déjà corrigée ; les
    cartes sont COMPOSÉES → l'image finale n'est rééchantillonnée qu'une fois.
    Retourne (map_x, map_y, gdx, gdy, n_valid, n_aps)."""
    H, W = ref.shape
    gdx, gdy, _ = global_shift(ref, mov, hann)
    dxs, dys, valid = ap_shifts(ref, mov, aps, box, gdx, gdy, max_local)
    map_x, map_y = build_maps((H, W), aps, dxs, dys, valid, gdx, gdy, grid_sigma)
    n_valid = int(valid.sum())
    for _ in range(max(0, n_pass - 1)):
        mov_w = cv2.remap(mov, map_x, map_y, cv2.INTER_CUBIC,
                          borderMode=cv2.BORDER_REPLICATE)
        rdx, rdy, rvalid = ap_shifts(ref, mov_w, aps, box, 0.0, 0.0, max_local)
        if int(rvalid.sum()) == 0:
            break
        rmap_x, rmap_y = build_maps((H, W), aps, rdx, rdy, rvalid, 0.0, 0.0, grid_sigma)
        map_x = cv2.remap(map_x, rmap_x, rmap_y, cv2.INTER_LINEAR,
                          borderMode=cv2.BORDER_REPLICATE)
        map_y = cv2.remap(map_y, rmap_x, rmap_y, cv2.INTER_LINEAR,
                          borderMode=cv2.BORDER_REPLICATE)
    return map_x, map_y, gdx, gdy, n_valid, len(aps)


# ══════════════════════════════════════════════════════════════════════════════
#  CŒUR ÉCLIPSE — détection du limbe lunaire (cercle RANSAC + Kåsa)
# ══════════════════════════════════════════════════════════════════════════════
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


# ── Siril helpers ─────────────────────────────────────────────────────────────
def get_siril_sequence(siril):
    """Détecte la séquence courante de Siril → dict ou None."""
    try:
        seq = siril.get_seq()
    except Exception:
        return None
    if seq is None:
        return None
    seqname = getattr(seq, "seqname", None)
    if not seqname:
        return None
    wd = os.getcwd()
    try:
        candidates = [f for f in Path(wd).iterdir()
                      if f.is_file() and f.suffix.lower() in FITS_EXT
                      and f.stem.startswith(seqname)]
    except Exception:
        candidates = []
    candidates.sort(key=lambda p: natural_key(p.name))
    files = [str(p) for p in candidates]
    number = int(getattr(seq, "number", len(files)) or 0)
    single_file = number > 1 and len(files) <= 1
    return {"seqname": seqname, "work_dir": wd, "files": files,
            "n_total": number, "nb_layers": int(getattr(seq, "nb_layers", 1) or 1),
            "single_file": single_file}

def siril_safe_cmd(siril, *args):
    try:
        siril.cmd(*args)
        return True
    except Exception:
        return False

def siril_safe_log(siril, msg, color=None):
    try:
        siril.log(msg, color=color) if color is not None else siril.log(msg)
    except Exception:
        pass


# ══════════════════════════════════════════════════════════════════════════════
#  Workers
# ══════════════════════════════════════════════════════════════════════════════
class RefLoadWorker(QThread):
    """Soleil : choisit et charge l'image de référence de la séquence."""
    progress = pyqtSignal(int, str)
    done     = pyqtSignal(object, object, int)   # lum, stretched, idx
    error    = pyqtSignal(str)

    def __init__(self, files, mode, manual_idx, parent=None):
        super().__init__(parent)
        self.files = files
        self.mode = mode
        self.manual_idx = manual_idx
        self._abort = False

    def abort(self): self._abort = True

    def run(self):
        files = self.files
        if not files:
            self.error.emit("Aucun fichier dans la séquence.")
            return
        n = len(files)
        if self.mode == "first":
            idx = 0
        elif self.mode == "last":
            idx = n - 1
        elif self.mode == "manual":
            idx = max(0, min(self.manual_idx, n - 1))
        else:   # auto — meilleure netteté
            best, best_idx = -1.0, 0
            for i, path in enumerate(files):
                if self._abort:
                    return
                try:
                    data, _ = load_fits(path)
                    score = estimate_sharpness(to_lum(data))
                    if score > best:
                        best, best_idx = score, i
                except Exception:
                    pass
                self.progress.emit(int(100*(i+1)/n),
                                   f"Analyse {i+1}/{n}  —  {Path(path).name}")
            idx = best_idx
        path = files[idx]
        self.progress.emit(99, f"Chargement référence : {Path(path).name}")
        try:
            data, _ = load_fits(path)
            self.done.emit(reg_lum(data), autostretch_lum(to_lum(data)), idx)
        except Exception as e:
            self.error.emit(f"Erreur chargement : {e}")


class MPAlignWorker(QThread):
    """Soleil : aligne toute la séquence (global + multi-points)."""
    progress = pyqtSignal(int, str)
    done     = pyqtSignal(int, int, str, bool)   # n_ok, n_total, out_dir, aborted
    error    = pyqtSignal(str)

    def __init__(self, params, parent=None):
        super().__init__(parent)
        self.params = params
        self._abort = False

    def abort(self): self._abort = True

    def _prepare_out_dir(self, out_folder, base):
        """Crée le sous-dossier isolé et purge nos sorties d'un run précédent."""
        try:
            d = Path(out_folder)
            d.mkdir(parents=True, exist_ok=True)
            for pat in (f"{base}_*.fit*", f"{base}.fit*", f"{base}.seq"):
                for old in d.glob(pat):
                    try: old.unlink()
                    except Exception: pass
            return True
        except Exception as e:
            self.error.emit(f"Création dossier sortie impossible : {e}")
            return False

    def run(self):
        p = self.params
        files = p["files"]; out_folder = p["out_folder"]; base = p["base"]
        ref = p["ref_lum"]; box = int(p["box"]); step = int(p["step"])
        max_local = float(p["max_local"]); n_pass = int(p["n_pass"])
        local_on = bool(p["local_on"])

        H, W = ref.shape
        n = len(files); n_ok = 0
        if not self._prepare_out_dir(out_folder, base):
            return

        def dst_for(i):
            return str(Path(out_folder) / f"{base}_{i+1:05d}.fit")

        hann = _hann((W, H))
        aps = []
        if local_on:
            self.progress.emit(0, "Points d'ancrage sur la référence…")
            aps = make_aps(ref, box, step)
            if aps:
                self.progress.emit(0, f"{len(aps)} points d'ancrage retenus "
                                      f"(boîte {box}px, pas {step}px).")
            else:
                self.progress.emit(0, "⚠ Aucun point d'ancrage (image trop lisse ?) "
                                      "→ translation globale seule.")

        for i, path in enumerate(files):
            if self._abort:
                break
            fname = Path(path).name
            try:
                data, hdr = load_fits(path)
                lum = reg_lum(data)
                if lum.shape != (H, W):
                    raise ValueError("dimensions différentes de la référence")
                if local_on and aps:
                    map_x, map_y, gdx, gdy, nv, na = mp_maps(
                        ref, lum, aps, box, max_local, hann, n_pass)
                    if abs(gdx) > W / 3 or abs(gdy) > H / 3:
                        raise ValueError(f"décalage global aberrant ({gdx:+.0f},{gdy:+.0f})")
                    aligned = warp_channels(data, map_x, map_y)
                    detail = f"global ({gdx:+.1f},{gdy:+.1f}) ; {nv}/{na} pts"
                else:
                    gdx, gdy, _ = global_shift(ref, lum, hann)
                    if abs(gdx) > W / 3 or abs(gdy) > H / 3:
                        raise ValueError(f"décalage global aberrant ({gdx:+.0f},{gdy:+.0f})")
                    M = np.float32([[1, 0, gdx], [0, 1, gdy]])
                    fl = cv2.INTER_CUBIC | cv2.WARP_INVERSE_MAP
                    if data.ndim == 2:
                        aligned = cv2.warpAffine(data, M, (W, H), flags=fl,
                                                 borderMode=cv2.BORDER_REPLICATE)
                    else:
                        aligned = np.stack([
                            cv2.warpAffine(data[c], M, (W, H), flags=fl,
                                           borderMode=cv2.BORDER_REPLICATE)
                            for c in range(data.shape[0])])
                    detail = f"global ({gdx:+.1f},{gdy:+.1f})"
                save_fits(aligned, hdr, dst_for(i))
                n_ok += 1
                self.progress.emit(int(100*(i+1)/n), f"✓ {i+1}/{n}  —  {fname}  [{detail}]")
            except Exception as e:
                self.progress.emit(int(100*(i+1)/n), f"✗ {i+1}/{n}  —  {fname} : {e}")
        self.done.emit(n_ok, n, out_folder, self._abort)


class DetectWorker(QThread):
    """Éclipse : détecte la Lune sur chaque pose."""
    result = pyqtSignal(int, object)   # index, (cx,cy,r,conf) ou None
    done   = pyqtSignal()

    def __init__(self, paths, parent=None):
        super().__init__(parent)
        self.paths = paths
        self._abort = False

    def abort(self): self._abort = True

    def run(self):
        for i, p in enumerate(self.paths):
            if self._abort:
                return
            try:
                data, _ = load_fits(p)
                res = detect_moon(to_lum(data))
            except Exception:
                res = None
            self.result.emit(i, res)
        self.done.emit()


class MoonAlignWorker(QThread):
    """Éclipse : recale les poses par translation (centres → référence)."""
    progress = pyqtSignal(int, str)
    done     = pyqtSignal(int, int, str)
    error    = pyqtSignal(str)

    def __init__(self, items, ref_cx, ref_cy, out_dir, scale=1.0, parent=None):
        super().__init__(parent)
        self.items = items          # [{path, cx, cy}]
        self.ref_cx = ref_cx
        self.ref_cy = ref_cy
        self.out_dir = out_dir
        self.scale = scale if scale > 0 else 1.0   # diviseur commun → sortie [0,1]
        self._abort = False

    def abort(self): self._abort = True

    def run(self):
        try:
            Path(self.out_dir).mkdir(parents=True, exist_ok=True)
        except Exception as e:
            self.error.emit(f"Dossier de sortie : {e}")
            return
        n = len(self.items); n_ok = 0
        for i, it in enumerate(self.items):
            if self._abort:
                break
            name = Path(it["path"]).name
            try:
                data, hdr = load_fits(it["path"])
                dy = self.ref_cy - it["cy"]
                dx = self.ref_cx - it["cx"]
                if data.ndim == 2:
                    out = ndshift(data, (dy, dx), order=1, mode="constant", cval=0.0)
                else:
                    out = ndshift(data, (0, dy, dx), order=1, mode="constant", cval=0.0)
                out = np.clip(out / self.scale, 0.0, 1.0)
                dst = str(Path(self.out_dir) / f"aligned_{n_ok+1:05d}.fit")
                save_fits(out, hdr, dst)
                n_ok += 1
                self.progress.emit(int(100*(i+1)/n),
                                   f"✓ {i+1}/{n}  —  {name}  (Δ {dx:+.1f},{dy:+.1f})")
            except Exception as e:
                self.progress.emit(int(100*(i+1)/n), f"✗ {i+1}/{n}  —  {name} : {e}")
        self.done.emit(n_ok, n, self.out_dir)


# ── Vue unifiée : pan/zoom, et (mode Lune) clic-glissé = centre, Maj+molette = rayon ─
class AlignView(QGraphicsView):
    center_at    = pyqtSignal(float, float)
    radius_delta = pyqtSignal(int)

    def __init__(self, scene, parent=None):
        super().__init__(scene, parent)
        self.setTransformationAnchor(QGraphicsView.ViewportAnchor.AnchorUnderMouse)
        self.setResizeAnchor(QGraphicsView.ViewportAnchor.AnchorUnderMouse)
        self.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform)
        self._pan = None
        self._dragging = False
        self._pick = False
        self.set_pick(False)

    def set_pick(self, on):
        """on=True (Lune) : clic-glissé pose le centre. on=False (Soleil) :
        clic-glissé = panoramique (ScrollHandDrag natif)."""
        self._pick = bool(on)
        self.setDragMode(QGraphicsView.DragMode.NoDrag if on
                         else QGraphicsView.DragMode.ScrollHandDrag)

    def wheelEvent(self, event):
        up = event.angleDelta().y() > 0
        if self._pick and (event.modifiers() & Qt.KeyboardModifier.ShiftModifier):
            self.radius_delta.emit(1 if up else -1)
        else:
            f = 1.15 if up else 1 / 1.15
            self.scale(f, f)

    def mousePressEvent(self, event):
        if self._pick and event.button() == Qt.MouseButton.LeftButton:
            self._dragging = True
            p = self.mapToScene(event.position().toPoint())
            self.center_at.emit(p.x(), p.y()); event.accept()
        elif event.button() == Qt.MouseButton.MiddleButton:
            self._pan = event.position().toPoint()
            self.setCursor(Qt.CursorShape.ClosedHandCursor); event.accept()
        else:
            super().mousePressEvent(event)

    def mouseMoveEvent(self, event):
        if self._dragging:
            p = self.mapToScene(event.position().toPoint())
            self.center_at.emit(p.x(), p.y()); event.accept()
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
#  Fenêtre principale — deux modes
# ══════════════════════════════════════════════════════════════════════════════
class AlignWindow(QMainWindow):

    def __init__(self, siril):
        super().__init__()
        self.siril = siril
        self.setWindowTitle(f"SirilJ Align — Soleil & Éclipse — v{VERSION}")
        self.setMinimumSize(1080, 690)
        self.setStyleSheet(DARK_SS)

        self.mode = "sun"
        self._last_dir = os.getcwd()
        self._H = self._W = 0
        self._pixitem = None

        # Soleil
        self.seq_info = None
        self.ref_lum = None; self.ref_stretched = None
        self.ref_idx = 0; self.ref_fname = ""
        self._out_subdir = ""; self._out_base = "mp"; self._aps = []
        self._ref_worker = None; self._align_worker = None
        self._ap_timer = QTimer(self); self._ap_timer.setSingleShot(True)
        self._ap_timer.setInterval(350); self._ap_timer.timeout.connect(self._refresh_aps)

        # Éclipse
        self._items = []; self._radius = 0.0; self._cur = -1; self._ref_idx = -1
        self._overlay = []; self._out_dir = ""
        self._detect = None; self._moon_align = None

        self._build_ui()
        self._refresh_sequence()

    # =========================================================================
    #  Construction de l'interface
    # =========================================================================
    def _build_ui(self):
        central = QWidget(); self.setCentralWidget(central)
        h = QHBoxLayout(central); h.setContentsMargins(6, 6, 6, 6); h.setSpacing(8)

        left = QWidget(); left.setFixedWidth(340)
        lv = QVBoxLayout(left); lv.setContentsMargins(0, 0, 0, 0); lv.setSpacing(6)

        trow = QHBoxLayout()
        trow.addWidget(QLabel("Cible :"))
        self.cmb_target = QComboBox()
        self.cmb_target.addItems(["☼ Soleil — surface (multi-points)",
                                  "☾ Éclipse — Lune (inter-poses)"])
        self.cmb_target.currentIndexChanged.connect(self._on_mode_changed)
        trow.addWidget(self.cmb_target, 1)
        lv.addLayout(trow)

        self.stack = QStackedWidget()
        self.stack.addWidget(self._build_sun_page())
        self.stack.addWidget(self._build_moon_page())
        lv.addWidget(self.stack, 1)
        h.addWidget(left)

        right = QWidget()
        rv = QVBoxLayout(right); rv.setContentsMargins(0, 0, 0, 0); rv.setSpacing(4)
        grp_prev = QGroupBox("Aperçu")
        pv = QVBoxLayout(grp_prev)
        self._scene = QGraphicsScene()
        self._view = AlignView(self._scene)
        self._view.setMinimumHeight(360)
        self._view.center_at.connect(self._on_center_at)
        self._view.radius_delta.connect(self._on_radius_delta)
        self.lbl_ref_info = QLabel("—")
        self.lbl_ref_info.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.lbl_ref_info.setStyleSheet("color:#888; font-size:10px;")
        pv.addWidget(self._view); pv.addWidget(self.lbl_ref_info)
        rv.addWidget(grp_prev, stretch=3)
        self.prog = QProgressBar(); rv.addWidget(self.prog)
        self.log = QTextEdit(); self.log.setReadOnly(True); self.log.setFixedHeight(110)
        rv.addWidget(self.log)
        h.addWidget(right, stretch=1)

        self.statusBar().showMessage("Prêt.")

    def _build_sun_page(self):
        page = QWidget(); lv = QVBoxLayout(page)
        lv.setContentsMargins(0, 0, 0, 0); lv.setSpacing(6)

        grp_seq = QGroupBox("Séquence Siril courante")
        sv = QVBoxLayout(grp_seq)
        self.lbl_seq_name = QLabel("—"); self.lbl_seq_name.setObjectName("SeqInfoKO")
        self.lbl_seq_name.setWordWrap(True)
        self.lbl_seq_info = QLabel(""); self.lbl_seq_info.setWordWrap(True)
        self.lbl_seq_info.setStyleSheet("color:#999; font-size:10px;")
        sv.addWidget(self.lbl_seq_name); sv.addWidget(self.lbl_seq_info)
        b_refresh = QPushButton("⟳  Rafraîchir"); b_refresh.clicked.connect(self._refresh_sequence)
        sv.addWidget(b_refresh)
        lv.addWidget(grp_seq)

        grp_ref = QGroupBox("Image de référence")
        rg = QGridLayout(grp_ref)
        rg.addWidget(QLabel("Mode :"), 0, 0)
        self.cmb_ref_mode = QComboBox()
        self.cmb_ref_mode.addItems(["Meilleure netteté", "Première", "Dernière", "Manuel"])
        self.cmb_ref_mode.currentIndexChanged.connect(
            lambda i: self.spn_ref_idx.setEnabled(i == 3))
        rg.addWidget(self.cmb_ref_mode, 0, 1)
        rg.addWidget(QLabel("Index :"), 1, 0)
        self.spn_ref_idx = QSpinBox(); self.spn_ref_idx.setRange(0, 99999)
        self.spn_ref_idx.setEnabled(False)
        rg.addWidget(self.spn_ref_idx, 1, 1)
        lv.addWidget(grp_ref)

        grp_mp = QGroupBox("Alignement multi-points")
        mg = QGridLayout(grp_mp)
        self.chk_local = QCheckBox("Correction locale (multi-points)")
        self.chk_local.setChecked(True)
        self.chk_local.setToolTip("Décoché : translation globale seule.")
        self.chk_local.toggled.connect(lambda _: self._schedule_aps())
        mg.addWidget(self.chk_local, 0, 0, 1, 2)
        mg.addWidget(QLabel("Boîte (px) :"), 1, 0)
        self.spn_box = QSpinBox(); self.spn_box.setRange(48, 384)
        self.spn_box.setSingleStep(16); self.spn_box.setValue(128)
        self.spn_box.setToolTip("Taille du patch de mesure. Petit = détail fin mais "
                                "fragile ; grand = robuste mais lisse.")
        self.spn_box.valueChanged.connect(lambda _: self._schedule_aps())
        mg.addWidget(self.spn_box, 1, 1)
        mg.addWidget(QLabel("Pas de grille (px) :"), 2, 0)
        self.spn_step = QSpinBox(); self.spn_step.setRange(16, 256)
        self.spn_step.setSingleStep(8); self.spn_step.setValue(64)
        self.spn_step.valueChanged.connect(lambda _: self._schedule_aps())
        mg.addWidget(self.spn_step, 2, 1)
        mg.addWidget(QLabel("Décalage local max (px) :"), 3, 0)
        self.spn_maxloc = QSpinBox(); self.spn_maxloc.setRange(2, 64); self.spn_maxloc.setValue(12)
        self.spn_maxloc.setToolTip("Au-delà, le point est rejeté (anti-aberration).")
        mg.addWidget(self.spn_maxloc, 3, 1)
        mg.addWidget(QLabel("Passes de raffinement :"), 4, 0)
        self.spn_pass = QSpinBox(); self.spn_pass.setRange(1, 4); self.spn_pass.setValue(2)
        mg.addWidget(self.spn_pass, 4, 1)
        self.lbl_aps = QLabel("points : —")
        self.lbl_aps.setStyleSheet("color:#999; font-size:10px;")
        mg.addWidget(self.lbl_aps, 5, 0, 1, 2)
        lv.addWidget(grp_mp)

        grp_out = QGroupBox("Sortie & intégration Siril")
        og = QGridLayout(grp_out)
        og.addWidget(QLabel("Préfixe :"), 0, 0)
        self.ed_prefix = QLineEdit("mp_")
        og.addWidget(self.ed_prefix, 0, 1)
        self.chk_convert = QCheckBox("Convertir en séquence Siril"); self.chk_convert.setChecked(True)
        og.addWidget(self.chk_convert, 1, 0, 1, 2)
        self.chk_stack = QCheckBox("Empiler après alignement"); self.chk_stack.setChecked(True)
        self.chk_stack.toggled.connect(lambda c: c and self.chk_convert.setChecked(True))
        og.addWidget(self.chk_stack, 2, 0, 1, 2)
        og.addWidget(QLabel("Méthode :"), 3, 0)
        self.cmb_stack = QComboBox(); self.cmb_stack.addItems(["sum", "med", "rej 3 3", "max", "min"])
        og.addWidget(self.cmb_stack, 3, 1)
        self.chk_load = QCheckBox("Charger le résultat dans Siril"); self.chk_load.setChecked(True)
        og.addWidget(self.chk_load, 4, 0, 1, 2)
        lv.addWidget(grp_out)

        hint = QLabel("Surface : clic-glissé = panoramique · molette = zoom. Les points "
                      "d'ancrage (bleus) sont posés sur la structure de la référence.")
        hint.setWordWrap(True); hint.setStyleSheet("font-size:8pt; color:#888;")
        lv.addWidget(hint)

        lv.addStretch()
        self.btn_prep = QPushButton("Préparer la référence"); self.btn_prep.setObjectName("BtnPrep")
        self.btn_prep.clicked.connect(self._prepare_ref)
        lv.addWidget(self.btn_prep)
        self.btn_start = QPushButton("▶  Démarrer l'alignement"); self.btn_start.setObjectName("BtnStart")
        self.btn_start.setEnabled(False); self.btn_start.clicked.connect(self._start_align)
        lv.addWidget(self.btn_start)
        self.btn_abort = QPushButton("⏹  Arrêter"); self.btn_abort.setEnabled(False)
        self.btn_abort.clicked.connect(self._abort_all)
        lv.addWidget(self.btn_abort)
        return page

    def _build_moon_page(self):
        page = QWidget(); lv = QVBoxLayout(page)
        lv.setContentsMargins(0, 0, 0, 0); lv.setSpacing(6)

        grp = QGroupBox("Poses (brackets ou masters par niveau)")
        gv = QVBoxLayout(grp)
        self.table = QTableWidget(0, 4)
        self.table.setHorizontalHeaderLabels(["Fichier", "X", "Y", "Conf."])
        self.table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        for c in (1, 2, 3):
            self.table.horizontalHeader().setSectionResizeMode(c, QHeaderView.ResizeMode.ResizeToContents)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setMaximumHeight(200)
        self.table.itemSelectionChanged.connect(self._on_row_selected)
        gv.addWidget(self.table)
        b_cur = QPushButton("⟳ Séquence Siril courante")
        b_cur.setToolTip("Charge les images de la séquence actuellement ouverte dans Siril "
                         "(fait automatiquement en entrant dans ce mode).")
        b_cur.clicked.connect(lambda: self._load_current_sequence(announce=True))
        gv.addWidget(b_cur)
        add_row = QHBoxLayout(); add_row.addWidget(QLabel("Ajouter :"))
        for txt, fn in (("FITS", self._add_files), ("Séquence", self._add_sequence),
                        ("Dossier", self._add_folder)):
            b = QPushButton(txt); b.clicked.connect(fn); add_row.addWidget(b)
        add_row.addStretch(); gv.addLayout(add_row)
        rm_row = QHBoxLayout(); rm_row.addWidget(QLabel("Retirer :"))
        for txt, fn in (("La sélection", self._remove_sel), ("Tout", self._clear)):
            b = QPushButton(txt); b.clicked.connect(fn); rm_row.addWidget(b)
        rm_row.addStretch(); gv.addLayout(rm_row)
        b_det = QPushButton("↻ Détecter la Lune (toutes)"); b_det.clicked.connect(self._detect_all)
        gv.addWidget(b_det)
        lv.addWidget(grp)

        grp2 = QGroupBox("Référence & rayon")
        g2 = QGridLayout(grp2)
        g2.addWidget(QLabel("Rayon lunaire (px) :"), 0, 0)
        self.spn_r = QSpinBox(); self.spn_r.setRange(0, 20000)
        self.spn_r.valueChanged.connect(self._on_radius_spin)
        g2.addWidget(self.spn_r, 0, 1)
        b_lock = QPushButton("Verrouiller au rayon médian"); b_lock.clicked.connect(self._lock_median_radius)
        g2.addWidget(b_lock, 1, 0, 1, 2)
        b_all = QPushButton("⊕ Appliquer le centre affiché à toutes les poses")
        b_all.setToolTip("La Lune bouge à peine dans un bracket : réglez un centre "
                         "parfait sur une pose nette, puis propagez-le.")
        b_all.clicked.connect(self._apply_center_all)
        g2.addWidget(b_all, 2, 0, 1, 2)
        self.lbl_ref = QLabel("Référence : meilleure confiance (auto)")
        self.lbl_ref.setStyleSheet("font-size:9pt; color:#aaa;"); self.lbl_ref.setWordWrap(True)
        g2.addWidget(self.lbl_ref, 3, 0, 1, 2)
        b_ref = QPushButton("Définir la pose affichée comme référence")
        b_ref.clicked.connect(self._set_ref_current)
        g2.addWidget(b_ref, 4, 0, 1, 2)
        lv.addWidget(grp2)

        grp3 = QGroupBox("Sortie & intégration Siril")
        g3 = QVBoxLayout(grp3)
        out_row = QHBoxLayout(); out_row.addWidget(QLabel("Dossier :"))
        self.ed_out = QLineEdit(); self.ed_out.setPlaceholderText("vide = dossier des poses")
        out_row.addWidget(self.ed_out)
        b_o = QPushButton("…"); b_o.setFixedWidth(28); b_o.clicked.connect(self._browse_out)
        out_row.addWidget(b_o); g3.addLayout(out_row)
        seq_row = QHBoxLayout(); seq_row.addWidget(QLabel("Séquence :"))
        self.chk_seq = QCheckBox("Créer aligned_.seq"); self.chk_seq.setChecked(True)
        seq_row.addWidget(self.chk_seq); seq_row.addStretch(); g3.addLayout(seq_row)
        lv.addWidget(grp3)

        hint = QLabel("Éclipse : clic-glissé = déplacer le centre (validé au relâchement) · "
                      "Maj+molette = rayon · molette = zoom · clic milieu = déplacer la vue. "
                      "Vérifiez le cercle bleu, surtout sur le diamant.")
        hint.setWordWrap(True); hint.setStyleSheet("font-size:8pt; color:#888;")
        lv.addWidget(hint)

        lv.addStretch()
        self.btn_align = QPushButton("▶  Aligner sur la Lune"); self.btn_align.setObjectName("BtnGo")
        self.btn_align.clicked.connect(self._do_align)
        lv.addWidget(self.btn_align)
        return page

    # =========================================================================
    #  Bascule de mode
    # =========================================================================
    def _on_mode_changed(self, idx):
        self.mode = "sun" if idx == 0 else "moon"
        self.stack.setCurrentIndex(idx)
        self._view.set_pick(self.mode == "moon")
        # réinitialise l'aperçu (les conventions d'affichage diffèrent)
        self._scene.clear(); self._pixitem = None; self._overlay = []
        self._H = self._W = 0
        if self.mode == "sun":
            self.lbl_ref_info.setText("Aucune référence chargée.")
            self._refresh_sequence()
        else:
            self.lbl_ref_info.setText("Sélectionnez une pose dans la liste.")
            if not self._items:                     # défaut : la séquence Siril courante
                self._load_current_sequence(announce=False)
            self._refresh_table()
        self.statusBar().showMessage("Mode Soleil." if self.mode == "sun" else "Mode Éclipse.")

    # =========================================================================
    #  ☼  SOLEIL
    # =========================================================================
    def _refresh_sequence(self):
        info = get_siril_sequence(self.siril)
        if info and info.get("single_file"):
            self.seq_info = None
            self.lbl_seq_name.setObjectName("SeqInfoKO")
            self.lbl_seq_name.setText(f"⚠ {info['seqname']} — séquence mono-fichier")
            self.lbl_seq_info.setText(
                f"FITSEQ/SER détecté ({info['n_total']} frames). Ce mode aligne des "
                "FITS individuels.\nExportez la séquence en FITS individuels puis ⟳.")
            self.btn_prep.setEnabled(False); self.btn_start.setEnabled(False)
        elif not info or not info["files"]:
            self.seq_info = None
            self.lbl_seq_name.setObjectName("SeqInfoKO")
            self.lbl_seq_name.setText("Aucune séquence chargée")
            self.lbl_seq_info.setText("Chargez une séquence dans Siril puis ⟳ Rafraîchir.")
            self.btn_prep.setEnabled(False); self.btn_start.setEnabled(False)
        else:
            self.seq_info = info
            self.lbl_seq_name.setObjectName("SeqInfoOK")
            self.lbl_seq_name.setText(f"✓  {info['seqname']}")
            layers = "monochrome" if info["nb_layers"] == 1 else "couleur (3 canaux)"
            self.lbl_seq_info.setText(f"{len(info['files'])} image(s)  •  {layers}\n"
                                      f"📂 {info['work_dir']}")
            self.spn_ref_idx.setRange(0, len(info["files"]) - 1)
            self.btn_prep.setEnabled(True)
            self._log(f"Séquence détectée : {info['seqname']} ({len(info['files'])} images)")
        self.lbl_seq_name.style().unpolish(self.lbl_seq_name)
        self.lbl_seq_name.style().polish(self.lbl_seq_name)

    def _prepare_ref(self):
        if not self.seq_info:
            self._log("⚠ Pas de séquence chargée.")
            return
        mode = ["auto", "first", "last", "manual"][self.cmb_ref_mode.currentIndex()]
        self.btn_prep.setEnabled(False); self.btn_start.setEnabled(False)
        self.btn_abort.setEnabled(True); self.prog.setValue(0)
        self._log("Recherche de la meilleure référence…" if mode == "auto"
                  else "Chargement de la référence…")
        self._ref_worker = RefLoadWorker(self.seq_info["files"], mode, self.spn_ref_idx.value())
        self._ref_worker.progress.connect(lambda p, m: (self.prog.setValue(p), self._log(m)))
        self._ref_worker.done.connect(self._on_ref_done)
        self._ref_worker.error.connect(lambda e: self._log(f"✗ {e}"))
        self._ref_worker.finished.connect(self._set_idle)
        self._ref_worker.start()

    def _on_ref_done(self, lum, stretched, idx):
        self.ref_lum = lum; self.ref_stretched = stretched; self.ref_idx = idx
        self.ref_fname = Path(self.seq_info["files"][idx]).name
        self.lbl_ref_info.setText(
            f"Référence #{idx}  —  {self.ref_fname}  |  {lum.shape[1]}×{lum.shape[0]} px")
        self.prog.setValue(100)
        self._log(f"✓ Référence prête : {self.ref_fname}")
        self.btn_start.setEnabled(True)
        self._refresh_aps()

    def _schedule_aps(self):
        if self.mode == "sun" and self.ref_lum is not None:
            self._ap_timer.start()

    def _refresh_aps(self):
        if self.ref_lum is None:
            return
        if self.chk_local.isChecked():
            self._aps = make_aps(self.ref_lum, self.spn_box.value(), self.spn_step.value())
            self.lbl_aps.setText(f"points : {len(self._aps)}")
        else:
            self._aps = []
            self.lbl_aps.setText("points : — (global seul)")
        self._draw_preview()

    def _draw_preview(self):
        """Référence + points d'ancrage gravés dans le pixmap (des milliers
        d'items de scène rendraient l'interface poussive)."""
        if self.ref_stretched is None:
            return
        self._H, self._W = self.ref_stretched.shape
        pm = lum_to_pixmap(self.ref_stretched)
        if self._aps:
            H = self.ref_stretched.shape[0]
            painter = QPainter(pm)
            painter.setPen(QPen(QColor("#00c8ff"), 2))
            for cy, cx in self._aps:
                sy = H - 1 - cy
                painter.drawRect(cx - 3, sy - 3, 6, 6)
            painter.end()
        if self._pixitem is None:
            self._pixitem = self._scene.addPixmap(pm)
        else:
            self._pixitem.setPixmap(pm)
        self._scene.setSceneRect(0, 0, pm.width(), pm.height())
        QTimer.singleShot(50, lambda: self._view.fitInView(
            self._scene.sceneRect(), Qt.AspectRatioMode.KeepAspectRatio))

    def _start_align(self):
        if not self.seq_info or self.ref_lum is None:
            return
        base = (self.ed_prefix.text().strip().rstrip("_") or "mp")
        out_subdir = str(Path(self.seq_info["work_dir"]) / base)
        self._out_subdir = out_subdir; self._out_base = base
        params = {"files": self.seq_info["files"], "out_folder": out_subdir, "base": base,
                  "ref_lum": self.ref_lum, "box": self.spn_box.value(),
                  "step": self.spn_step.value(), "max_local": self.spn_maxloc.value(),
                  "n_pass": self.spn_pass.value(), "local_on": self.chk_local.isChecked()}
        self.btn_start.setEnabled(False); self.btn_prep.setEnabled(False)
        self.btn_abort.setEnabled(True); self.prog.setValue(0)
        mode_txt = "multi-points" if params["local_on"] else "global seul"
        self._log(f"━━ Alignement {mode_txt} → {out_subdir}")
        siril_safe_log(self.siril, f"SirilJ Align — alignement {mode_txt} (séquence : {base})")
        self._align_worker = MPAlignWorker(params)
        self._align_worker.progress.connect(lambda p, m: (self.prog.setValue(p), self._log(m)))
        self._align_worker.done.connect(self._on_sun_done)
        self._align_worker.error.connect(lambda e: self._log(f"✗ {e}"))
        self._align_worker.finished.connect(self._set_idle)
        self._align_worker.start()

    def _on_sun_done(self, n_ok, n, out_dir, aborted):
        self.prog.setValue(100)
        if aborted:
            self._log(f"⏹ Interrompu : {n_ok}/{n} images écrites — "
                      "étapes Siril (convert/stack/load) non lancées.")
            return
        msg = f"✓ Alignement terminé : {n_ok}/{n} images."
        self._log(msg); siril_safe_log(self.siril, msg); self.statusBar().showMessage(msg)
        if n_ok == 0:
            self._log("⚠ Aucun fichier produit — étapes Siril non lancées.")
            return
        if not (self.chk_convert.isChecked() or self.chk_stack.isChecked()
                or self.chk_load.isChecked()):
            return
        base = self._out_base
        subdir_abs = str(Path(out_dir).resolve())
        orig_wd = os.getcwd()
        if not siril_safe_cmd(self.siril, f'cd "{subdir_abs}"'):   # guillemets : espaces
            self._log(f"⚠ Impossible de faire 'cd {subdir_abs}' dans Siril. Étapes annulées.")
            return
        try:
            if self.chk_convert.isChecked() or self.chk_stack.isChecked():
                self._log(f"→ Conversion en séquence Siril : {base}")
                ok = siril_safe_cmd(self.siril, "convert", base, "-fitseq")
                if not ok:
                    ok = siril_safe_cmd(self.siril, "convert", base)
                self._log(f"✓ Séquence «{base}» créée." if ok else "⚠ Échec de convert.")
            result_file = None
            if self.chk_stack.isChecked():
                method = self.cmb_stack.currentText()
                out_img = f"{base}_stacked"
                self._log(f"→ Empilement ({method}) → {out_img}")
                args = ["stack", base] + method.split() + [f"-out={out_img}"]
                if siril_safe_cmd(self.siril, *args):
                    self._log(f"✓ Empilement terminé : {out_img}"); result_file = out_img
                else:
                    self._log("⚠ Échec de la commande stack.")
            if self.chk_load.isChecked():
                if result_file is None:
                    result_file = f"{base}_00001"
                self._log(f"→ Chargement dans Siril : {result_file}")
                if siril_safe_cmd(self.siril, "load", result_file):
                    self._log("✓ Affichage dans Siril mis à jour.")
                else:
                    self._log(f"⚠ Échec de load (chargez {result_file} manuellement).")
        finally:
            siril_safe_cmd(self.siril, f'cd "{orig_wd}"')   # guillemets : espaces

    def _set_idle(self):
        self.btn_prep.setEnabled(self.seq_info is not None)
        self.btn_abort.setEnabled(False)
        self.btn_start.setEnabled(self.ref_lum is not None)

    def _abort_all(self, wait=False):
        for w in (self._ref_worker, self._align_worker, self._detect, self._moon_align):
            if w and w.isRunning():
                w.abort()
                if wait:
                    w.wait(5000)
        if not wait:
            self._log("⏹ Arrêté.")

    # =========================================================================
    #  ☾  ÉCLIPSE — liste des poses
    # =========================================================================
    def _load_current_sequence(self, announce=True):
        """Charge dans la table les images de la séquence ouverte dans Siril."""
        info = get_siril_sequence(self.siril)
        if info and info.get("single_file"):
            if announce:
                self._log(f"⚠ {info['seqname']} : séquence mono-fichier (FITSEQ/SER). "
                          "Exportez en FITS individuels, ou utilisez « Ajouter ».")
            return False
        if not info or not info["files"]:
            if announce:
                self._log("Aucune séquence Siril courante — utilisez « Ajouter ».")
            return False
        self._items.clear(); self._cur = -1; self._ref_idx = -1
        self._add_paths(info["files"])
        self._log(f"Séquence Siril «{info['seqname']}» chargée "
                  f"({len(info['files'])} poses). Lancez « Détecter la Lune ».")
        return True

    def _add_paths(self, paths):
        added = 0
        for p in paths:
            try:
                data, _ = load_fits(p)
                self._items.append({"path": p, "cx": data.shape[-1] / 2.0,
                                    "cy": data.shape[-2] / 2.0, "r": 0.0, "conf": 0.0,
                                    "shape": data.shape, "dmax": float(np.max(data))})
                added += 1
            except Exception as e:
                self._log(f"✗ {Path(p).name} : {e}")
        if added:
            self._refresh_table()
            self._log(f"+ {added} pose(s). Lancez « Détecter la Lune ».")

    def _add_files(self):
        files, _ = QFileDialog.getOpenFileNames(self, "Ajouter des poses", self._last_dir,
            "FITS (*.fit *.fits *.fts *.FIT *.FITS *.FTS)")
        if files:
            self._last_dir = str(Path(files[0]).parent)
            self._add_paths(files)

    def _add_folder(self):
        d = QFileDialog.getExistingDirectory(self, "Dossier", self._last_dir)
        if not d:
            return
        self._last_dir = d
        files = sorted(str(f) for f in Path(d).iterdir()
                       if f.is_file() and f.suffix.lower() in FITS_EXT)
        if files:
            self._add_paths(files)
        else:
            self._log("Aucun FITS dans ce dossier.")

    def _add_sequence(self):
        f, _ = QFileDialog.getOpenFileName(self, "Séquence Siril", self._last_dir,
                                           "Séquence Siril (*.seq)")
        if not f:
            return
        self._last_dir = str(Path(f).parent)
        seqname = None
        try:
            with open(f, encoding="utf-8", errors="ignore") as fh:
                for line in fh:
                    m = re.match(r"\s*S\s+'([^']*)'", line)
                    if m:
                        seqname = m.group(1); break
        except Exception as e:
            self._log(f"✗ Lecture .seq : {e}"); return
        if not seqname:
            self._log("Nom de séquence introuvable dans le .seq."); return
        folder = Path(f).parent
        files = sorted(str(p) for p in folder.iterdir()
                       if p.is_file() and p.suffix.lower() in FITS_EXT
                       and p.stem.startswith(seqname))
        if files:
            self._add_paths(files)
            self._log(f"Séquence «{seqname}» : {len(files)} image(s) chargée(s).")
        else:
            self._log(f"Aucun FITS «{seqname}*» (séquence mono-fichier FITSEQ/SER ?).")

    def _remove_sel(self):
        rows = sorted({i.row() for i in self.table.selectedIndexes()}, reverse=True)
        for r in rows:
            if 0 <= r < len(self._items):
                del self._items[r]
        self._cur = -1
        self._refresh_table()

    def _clear(self):
        self._items.clear(); self._cur = -1; self._ref_idx = -1
        self._refresh_table()

    def _refresh_table(self):
        self.table.blockSignals(True)
        self.table.setRowCount(len(self._items))
        for i, it in enumerate(self._items):
            nm = QTableWidgetItem(Path(it["path"]).name); nm.setToolTip(it["path"])
            cx = QTableWidgetItem(f"{it['cx']:.0f}")
            cy = QTableWidgetItem(f"{it['cy']:.0f}")
            cf = QTableWidgetItem(f"{it['conf']:.2f}" if it["conf"] else "—")
            if it["conf"]:
                col = QColor("#7ec77e") if it["conf"] >= 0.15 else \
                      (QColor("#d9b34a") if it["conf"] >= 0.07 else QColor("#d77"))
                cf.setForeground(col)
            for c, itm in enumerate((nm, cx, cy, cf)):
                itm.setFlags(itm.flags() & ~Qt.ItemFlag.ItemIsEditable)
                self.table.setItem(i, c, itm)
        self.table.blockSignals(False)

    # =========================================================================
    #  ☾  ÉCLIPSE — détection
    # =========================================================================
    def _detect_all(self):
        if not self._items:
            return
        self.prog.setValue(0)
        self._log("Détection du limbe (RANSAC)…")
        self._detect = DetectWorker([it["path"] for it in self._items])
        self._detect.result.connect(self._on_detect_one)
        self._detect.done.connect(self._on_detect_done)
        self._detect.start()

    def _on_detect_one(self, i, res):
        if 0 <= i < len(self._items) and res is not None:
            cx, cy, r, conf = res
            self._items[i].update(cx=cx, cy=cy, r=r, conf=conf)
        self.prog.setValue(int(100 * (i + 1) / max(len(self._items), 1)))

    def _on_detect_done(self):
        self._refresh_table()
        self._lock_median_radius()
        if self._ref_idx < 0:
            self._ref_idx = max(range(len(self._items)),
                                key=lambda k: self._items[k]["conf"], default=-1)
            self._update_ref_label()
        n_ok = sum(1 for it in self._items if it["conf"] > 0)
        self._log(f"✓ Détection : {n_ok}/{len(self._items)} poses. "
                  "Vérifiez les faibles confiances (orange/rouge) au clic.")
        if self._items and self._cur < 0:
            self._select_row(0)

    def _lock_median_radius(self):
        rs = [it["r"] for it in self._items if it["conf"] > 0 and it["r"] > 0]
        if rs:
            self._radius = float(np.median(rs))
            self.spn_r.blockSignals(True); self.spn_r.setValue(int(round(self._radius)))
            self.spn_r.blockSignals(False)
            self._draw_overlay()
            self._log(f"Rayon lunaire verrouillé à {self._radius:.0f} px (médiane).")

    # =========================================================================
    #  ☾  ÉCLIPSE — aperçu & édition manuelle
    # =========================================================================
    def _on_row_selected(self):
        rows = [i.row() for i in self.table.selectedIndexes()]
        if rows:
            self._select_row(rows[0])

    def _select_row(self, idx):
        if not (0 <= idx < len(self._items)):
            return
        self._cur = idx
        try:
            data, _ = load_fits(self._items[idx]["path"])
        except Exception as e:
            self._log(f"✗ {e}"); return
        lum = to_lum(data)
        self._H, self._W = lum.shape
        pm = lum_to_pixmap(autostretch_lum(lum))
        if self._pixitem is None:
            self._pixitem = self._scene.addPixmap(pm)
        else:
            self._pixitem.setPixmap(pm)
        self._scene.setSceneRect(0, 0, pm.width(), pm.height())
        QTimer.singleShot(30, lambda: self._view.fitInView(
            self._scene.sceneRect(), Qt.AspectRatioMode.KeepAspectRatio))
        self._draw_overlay()

    def _draw_overlay(self):
        for it in self._overlay:
            self._scene.removeItem(it)
        self._overlay.clear()
        if self.mode != "moon" or self._cur < 0 or self._H == 0:
            return
        it = self._items[self._cur]
        sx = it["cx"]; sy = self._H - 1 - it["cy"]
        pen = QPen(QColor("#ff3333")); pen.setCosmetic(True); pen.setWidth(2)
        L = 12
        self._overlay.append(self._scene.addLine(sx - L, sy, sx + L, sy, pen))
        self._overlay.append(self._scene.addLine(sx, sy - L, sx, sy + L, pen))
        lr = self._radius if self._radius > 0 else it["r"]
        if lr > 0:
            penc = QPen(QColor("#33aaff")); penc.setCosmetic(True); penc.setWidth(2)
            self._overlay.append(self._scene.addEllipse(sx - lr, sy - lr, 2*lr, 2*lr, penc))

    def _on_center_at(self, sx, sy):
        # mode Lune uniquement (pick actif) : le centre suit la souris.
        if self.mode != "moon" or self._cur < 0 or self._H == 0:
            return
        it = self._items[self._cur]
        it["cx"] = float(sx); it["cy"] = float(self._H - 1 - sy); it["conf"] = 1.0
        self._draw_overlay()
        self.table.blockSignals(True)
        for col, txt in ((1, f"{it['cx']:.0f}"), (2, f"{it['cy']:.0f}"), (3, "1.00")):
            cell = self.table.item(self._cur, col)
            if cell:
                cell.setText(txt)
                if col == 3:
                    cell.setForeground(QColor("#7ec77e"))
        self.table.blockSignals(False)

    def _on_radius_delta(self, direction):
        step = max(1, int(self.spn_r.value() * 0.03))
        self.spn_r.setValue(self.spn_r.value() + direction * step)

    def _on_radius_spin(self, val):
        self._radius = float(val)
        self._draw_overlay()

    # =========================================================================
    #  ☾  ÉCLIPSE — référence & alignement
    # =========================================================================
    def _apply_center_all(self):
        if self._cur < 0:
            self._log("Affichez d'abord une pose et réglez son centre.")
            return
        cx = self._items[self._cur]["cx"]; cy = self._items[self._cur]["cy"]
        for it in self._items:
            it["cx"] = cx; it["cy"] = cy; it["conf"] = 1.0
        self._refresh_table()
        self.table.selectRow(self._cur)
        self._draw_overlay()
        self._log(f"Centre ({cx:.0f}, {cy:.0f}) appliqué à toutes les poses.")

    def _set_ref_current(self):
        if self._cur >= 0:
            self._ref_idx = self._cur
            self._update_ref_label()

    def _update_ref_label(self):
        if 0 <= self._ref_idx < len(self._items):
            self.lbl_ref.setText(f"Référence : {Path(self._items[self._ref_idx]['path']).name}")

    def _do_align(self):
        if len(self._items) < 2:
            self._log("Ajoutez au moins deux poses."); return
        if self._ref_idx < 0:
            self._ref_idx = max(range(len(self._items)),
                                key=lambda k: self._items[k]["conf"], default=0)
            self._update_ref_label()
        ref = self._items[self._ref_idx]
        out = self.ed_out.text().strip() or str(Path(self._items[0]["path"]).parent)
        self._out_dir = out
        gmax = max((it.get("dmax", 1.0) for it in self._items), default=1.0)
        scale = gmax if gmax > 1.5 else 1.0
        self.btn_align.setEnabled(False); self.prog.setValue(0)
        self._log(f"━━ Alignement sur la Lune → {out}  "
                  f"(réf : {Path(ref['path']).name}, échelle ÷{scale:.0f})")
        self._moon_align = MoonAlignWorker(
            [{"path": it["path"], "cx": it["cx"], "cy": it["cy"]} for it in self._items],
            ref["cx"], ref["cy"], out, scale=scale)
        self._moon_align.progress.connect(lambda v, m: (self.prog.setValue(v), self._log(m)))
        self._moon_align.done.connect(self._on_moon_done)
        self._moon_align.error.connect(lambda e: (self._log(f"✗ {e}"),
                                                  self.btn_align.setEnabled(True)))
        self._moon_align.start()

    def _on_moon_done(self, n_ok, n, out_dir):
        self.btn_align.setEnabled(True)
        self._log(f"✓ {n_ok}/{n} poses alignées → {out_dir}")
        siril_safe_log(self.siril, f"SirilJ Align — {n_ok}/{n} poses alignées sur la Lune.")
        if self.chk_seq.isChecked() and n_ok > 0:
            self._write_seq(out_dir, n_ok)

    def _write_seq(self, out_dir, n_ok):
        """Écrit aligned_.seq directement (ne référence que aligned_NNNNN.fit,
        donc aucun mélange même si d'autres FITS traînent — contrairement à convert)."""
        layers = 3 if (self._items and len(self._items[0]["shape"]) == 3) else 1
        ref0 = max(0, min(self._ref_idx, n_ok - 1))
        try:
            lines = ["#Siril sequence file. Généré par SirilJ Align.",
                     f"S 'aligned_' 1 {n_ok} {n_ok} 5 {ref0} 1", f"L {layers}"]
            lines += [f"I {i} 1" for i in range(1, n_ok + 1)]
            (Path(out_dir) / "aligned_.seq").write_text("\n".join(lines) + "\n")
            self._log("✓ Séquence «aligned_.seq» créée. "
                      "Dans Siril : Ouvrir une séquence → aligned_.seq")
        except Exception as e:
            self._log(f"⚠ aligned_.seq : {e}")

    def _browse_out(self):
        d = QFileDialog.getExistingDirectory(self, "Dossier de sortie", self._last_dir)
        if d:
            self._last_dir = d
            self.ed_out.setText(d)

    # =========================================================================
    def _log(self, msg):
        esc = html.escape(msg)
        weight = "bold" if msg.lstrip().startswith(("⚠", "✗")) else "normal"
        self.log.append(f'<div style="font-weight:{weight}">{esc}</div>')
        self.log.ensureCursorVisible()
        self.statusBar().showMessage(msg)

    def closeEvent(self, event):
        self._abort_all(wait=True)
        super().closeEvent(event)


def main():
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    siril = s.SirilInterface()
    try:
        siril.connect()
    except Exception as e:
        print(f"Avertissement : connexion Siril impossible ({e}).")
    win = AlignWindow(siril)
    win.show()
    app.exec()


if __name__ == "__main__":
    main()
