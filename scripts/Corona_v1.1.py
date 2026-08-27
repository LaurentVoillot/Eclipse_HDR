##############################################
# Corona v1.1 — détail tangentiel (ACHF)
# Rehaussement de couronne solaire (éclipse)
# Version 1.1.0 — variante expérimentale (basée sur Corona 3.0.0)
##############################################
#
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Le script original Corona.py reste inchangé — ceci est une variante à
# part qui ajoute la méthode « Détail tangentiel » (ACHF simplifié).
#
# Crédits
# -------
#   • Détail tangentiel — masque flou le long des arcs de cercle, ACHF
#     simplifié (d'après Druckmüller, Contrib. Astron. Obs. Skalnaté
#     Pleso 2006 — méthode de ses images d'éclipse de référence)
#   • RHEF — Radial Histogram Equalizing Filter
#     (Gilly & Cranmer, Solar Physics 2025)
#   • FNRGF — Fourier Normalizing Radial-Graded Filter
#     (Druckmüllerová, Morgan & Habbal, ApJ 2011)
#   • MGN — Multi-Scale Gaussian Normalization
#     (Morgan & Druckmüller, Solar Physics 2014)
#   • Interface : PyQt6 / conventions VeraLux

"""
Corona — révèle la structure de la couronne solaire (streamers, jets,
boucles). Trois méthodes éprouvées (références sunkit-image / SunPy) :

RHEF — Radial Histogram Equalizing Filter (Gilly & Cranmer 2025) [défaut]
    Chaque pixel est remplacé par son rang en percentile DANS son anneau
    radial. Profil de sortie quasi uniforme, SANS réglage, sans fond gris
    ni lobes sombres. Idéal comme image à finir au GHS dans Siril.

FNRGF — Fourier Normalizing Radial-Graded Filter (Druckmüllerová 2011)
    (I − moyenne) / écart-type, où moyenne ET écart-type locaux sont
    approchés par une série de Fourier azimutale d'ordre n → l'« attendu »
    suit les streamers (pas de lobes sombres du NRGF simple). `ordre`
    faible = grandes structures (naturel) ; élevé = très aplati.

MGN — Multi-Scale Gaussian Normalization (Morgan & Druckmüller 2014)
    détail = arctan(k · (image − C_σ) / S_σ) moyenné sur plusieurs σ,
    mélangé à une tonalité globale (gamma) pondérée par h. Purement local,
    sans centre ni masque.

Détail tangentiel — ACHF simplifié (Druckmüller 2006) [nouveau en v1.1]
    Masque flou calculé LE LONG des arcs de cercle (jamais en radial),
    soustrait de l'image : seules les structures RADIALES (streamers,
    jets, plumes polaires) ressortent, le gradient radial reste intact.
    Exploite la géométrie connue de la couronne — ce que les filtres
    statistiques ne font pas. Se combine bien : RHEF d'abord (aplatir),
    recharger le résultat, puis Détail tangentiel (affûter).

Entrée : l'image courante de Siril (typiquement le HDR de FusionHDR) ou
un fichier FITS. Données supposées linéaires.

Compatibilité
-------------
• Siril 1.3+   • Python 3.10+ (via sirilpy)
• Dépendances : numpy, astropy, scipy, PyQt6
"""

import sys
import os
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
                           binary_closing, map_coordinates)

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QGridLayout, QLabel, QPushButton, QComboBox, QCheckBox, QDoubleSpinBox,
    QSpinBox, QProgressBar, QTextEdit, QFileDialog, QGroupBox,
    QGraphicsView, QGraphicsScene, QGraphicsPixmapItem
)
from PyQt6.QtCore import Qt, QThread, pyqtSignal, QTimer
from PyQt6.QtGui import QImage, QPixmap, QPainter, QColor, QPen

VERSION = "1.1.0"

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
QComboBox          { background:#3c3c3c; border:1px solid #555; border-radius:3px;
                     padding:3px 6px; color:#d4d4d4; }
QComboBox QAbstractItemView { background:#2b2b2b; color:#d4d4d4;
                               selection-background-color:#555; }
QDoubleSpinBox     { background:#1e1e1e; border:1px solid #555; border-radius:3px;
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
QScrollBar:vertical   { background:#2b2b2b; width:10px; }
QScrollBar::handle:vertical { background:#555; border-radius:5px; }
"""

# ── Helpers ──────────────────────────────────────────────────────────────────
def _mtf(x, m, lo, hi):
    dist = hi - lo
    if dist < 1e-9:
        return np.where(x > lo, 1.0, 0.0).astype(np.float32)
    xp = np.clip((x - lo) / dist, 0.0, 1.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        y = ((m - 1.0) * xp) / ((2.0 * m - 1.0) * xp - m)
    return np.nan_to_num(y, nan=0.0, posinf=1.0, neginf=0.0).astype(np.float32)

def to_channels(d):
    return [d] if d.ndim == 2 else [d[i] for i in range(d.shape[0])]

def autostretch_display(data):
    chans = to_channels(data)
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

def neutralize_corona(img, amount, balance=True):
    """Désature la couronne en préservant le rouge Hα des protubérances.

    Corona ne CRÉE pas de saturation — elle multiplie les trois canaux par le
    même facteur, donc les rapports R:V:B sont conservés. Mais en remontant la
    couronne faible à pleine luminosité, elle RÉVÈLE une dominante jusque-là
    invisible. Or la couronne K est de la lumière photosphérique diffusée par
    les électrons : elle est physiquement BLANCHE, sa teinte est un artefact.

    La neutraliser maximise le RAPPORT de saturation avec les protubérances —
    la saturation absolue, elle, se remonte ensuite d'un curseur.
    Aucun centre n'est requis (fonctionne donc aussi en mode MGN).
    amount : 0 = inchangé, 1 = couronne parfaitement neutre.
    """
    # La balance des blancs doit s'appliquer même si la désaturation est à 0 :
    # ce sont deux réglages indépendants (retirer la dominante ≠ désaturer).
    if img.ndim != 3 or img.shape[0] < 3 or (amount <= 1e-6 and not balance):
        return img.astype(np.float32), (1.0, 1.0, 1.0)
    out = img.astype(np.float32).copy()
    lum = out.mean(axis=0)
    # Zone de référence : la moitié la plus lumineuse (exclut disque et ciel
    # vide, dont les rapports de couleur sont dominés par le bruit).
    thr = float(np.percentile(lum, 50.0))
    zone = lum > max(thr, 1e-9)
    if int(zone.sum()) < 500:
        zone = lum > 0

    def _ratio(im):
        """R / max(V, B) — critère ABSOLU de « vraiment rouge ».

        Un seuil relatif (rougeur au-dessus de la médiane) protège tout ce qui
        est un peu plus chaud que la moyenne : une dominante jaune passe pour
        du Hα, garde sa saturation pendant que le reste est désaturé, et
        ressort donc en taches jaunes au lieu de disparaître.
        Le rapport R/max(V,B) tranche sans ambiguïté : une couronne neutre vaut
        ~1,0 · une dominante jaune ~1,1 · une protubérance Hα de 3 à 10."""
        return im[0] / np.maximum(np.maximum(im[1], im[2]), 1e-12)

    def _hw(im):
        t = np.clip((_ratio(im) - 1.15) / (1.35 - 1.15), 0.0, 1.0)
        return (t * t * (3.0 - 2.0 * t)).astype(np.float32)

    # Poids Hα calculé sur l'image D'ORIGINE, avant toute balance : celle-ci
    # remonte le canal le plus faible (souvent le bleu) pour neutraliser la
    # couronne, ce qui écrase le rapport R/max(V,B) des protubérances et leur
    # ferait perdre leur protection.
    hw = _hw(out)

    gains = (1.0, 1.0, 1.0)
    if balance:
        neutral = zone & (hw < 0.2)                  # exclut les protubérances
        ref = neutral if int(neutral.sum()) > 500 else zone
        meds = [float(np.median(out[c][ref])) for c in range(3)]
        tgt = float(np.mean(meds))
        if all(m > 1e-9 for m in meds):
            gains = tuple(tgt / m for m in meds)
            for c in range(3):
                out[c] = out[c] * gains[c]
    s = (1.0 - float(amount) * (1.0 - hw)).astype(np.float32)
    lum = out.mean(axis=0)
    out = np.maximum(lum[None, :, :] + (out - lum[None, :, :]) * s[None, :, :], 0.0)
    return out.astype(np.float32), gains


def save_fits(data, path):
    fits.PrimaryHDU(np.asarray(data, dtype=np.float32)).writeto(path, overwrite=True)


# ── MGN — Multi-Scale Gaussian Normalization ────────────────────────────────────
def mgn(lum, scales=(1.25, 2.5, 5, 10, 20, 40), k=0.7, gamma=3.2, h=0.7,
        denoise=2.0, progress=None):
    """Normalisation gaussienne multi-échelle (Morgan & Druckmüller 2014).

    lum : luminance 2-D. Retourne une image rehaussée ~[0,1].
    denoise : plancher de bruit (× σ bruit estimé) sur l'écart-type local —
    empêche d'amplifier le grain dans les zones plates/faibles. 0 = désactivé."""
    a = lum.astype(np.float32)
    amin, amax = float(a.min()), float(a.max())
    g = ((a - amin) / (amax - amin + 1e-9)) ** (1.0 / max(gamma, 1e-3))  # tonalité globale

    # Plancher = denoise × σ_bruit (MAD du résidu haute fréquence). Dans les
    # zones plates, std ≈ bruit ; sans plancher, MGN divise par le bruit et le
    # remonte à plein contraste (grain). Avec plancher, ces zones sont calmées
    # tandis que la vraie structure (std ≫ plancher) reste normalisée.
    floor = 0.0
    if denoise > 0:
        hi = a - gaussian_filter(a, 1.0, truncate=3.0)
        sigma_n = 1.4826 * float(np.median(np.abs(hi - np.median(hi))))
        floor = float(denoise) * sigma_n

    detail = np.zeros_like(a)
    n = len(scales)
    for i, sgm in enumerate(scales):
        mean = gaussian_filter(a, sgm, truncate=3.0)
        var = gaussian_filter((a - mean) ** 2, sgm, truncate=3.0)
        std = np.sqrt(np.maximum(var, 0.0))
        std = np.maximum(std, floor)              # plancher de bruit
        std = np.where(std < 1e-6, 1e-6, std)
        detail += np.arctan(k * (a - mean) / std)
        if progress:
            progress(20 + int(60 * (i + 1) / n), f"MGN — échelle {sgm:g} px")
    detail /= n
    # détail arctan ∈ [-π/2, π/2] ramené en ~[0,1], mélangé à la tonalité globale
    return (h * g + (1.0 - h) * (0.5 + detail / np.pi)).astype(np.float32)


# ── Détection du disque (centroïde du disque sombre enclos) ──────────────────────
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

def detect_disk(lum):
    """Centroïde du disque sombre enclos par la couronne → (cx, cy, r) ou None.
    Fiable ici car Corona traite le HDR fusionné (couronne brillante tout autour)."""
    a = gaussian_filter(lum.astype(np.float32), 2.0)
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


# ── NRGF — Normalized Radial Gradient Filter (Morgan 2006) ───────────────────────
def _noise_sigma(a):
    """Estimation robuste du bruit : 1,4826 × MAD(I − gaussien(I, 1))."""
    hi = a - gaussian_filter(a, 1.0, truncate=3.0)
    return 1.4826 * float(np.median(np.abs(hi - np.median(hi))))


def rhef(lum, cx, cy, lunar_r=0.0, denoise=2.0, rmax_px=0.0):
    """RHEF — Radial Histogram Equalizing Filter (Gilly & Cranmer, Solar Phys. 2025).
    Chaque pixel est remplacé par son rang en percentile DANS son anneau radial
    → profil de sortie quasi uniforme, sans réglage. Pas de fond gris ni de lobes
    sombres (transformation monotone, bornée [0,1]). Disque + au-delà → sombre.

    Garde-fou : l'égalisation d'histogramme étire n'importe quel signal à [0,1],
    donc un anneau de ciel vide deviendrait du bruit plein contraste. On pondère
    donc chaque anneau par sa confiance signal (écart-type / bruit) : les anneaux
    dominés par le bruit restent sombres, ceux avec vrai signal sont pleinement
    égalisés."""
    a = lum.astype(np.float32)
    sig0 = _noise_sigma(a)                # bruit par pixel AVANT lissage (référence)
    red = 1.0
    if denoise > 0:                       # léger lissage : limite l'amplification du grain
        ss = 0.5 + 0.4 * float(denoise)
        a = gaussian_filter(a, ss, truncate=3.0)
        red = max(2.0 * ss * np.sqrt(np.pi), 1.0)   # réduction du bruit blanc par le lissage
    H, W = a.shape
    yy, xx = np.indices((H, W)).astype(np.float32)
    r = np.hypot(xx - cx, yy - cy)
    rmax = float(rmax_px) if rmax_px and rmax_px > 0 else float(r.max())
    nb = int(r.max()) + 1
    ri = np.clip(r.astype(np.int32), 0, nb - 1)
    keep = (r >= lunar_r) & (r <= rmax) & np.isfinite(a)
    out = np.zeros((H, W), np.float32)
    kr = ri[keep]
    kv = a[keep].astype(np.float64)
    n = kr.size
    if n == 0:
        return out
    cnt = np.bincount(kr, minlength=nb)
    # rang dans l'anneau = position dans le tri (anneau, valeur) − début de l'anneau
    order = np.lexsort((kv, kr))          # tri primaire anneau, secondaire valeur
    gpos = np.empty(n, dtype=np.int64)
    gpos[order] = np.arange(n)
    ring_start = np.concatenate(([0], np.cumsum(cnt)[:-1]))
    rank = (gpos - ring_start[kr]) / np.maximum(cnt[kr] - 1, 1)   # [0,1] par anneau
    # poids de confiance signal par anneau (anti-bruit du bord vide)
    g = cnt > 0
    mr = np.zeros(nb); std_r = np.zeros(nb)
    mr[g] = np.bincount(kr, weights=kv, minlength=nb)[g] / cnt[g]
    std_r[g] = np.sqrt(np.maximum(
        np.bincount(kr, weights=kv * kv, minlength=nb)[g] / cnt[g] - mr[g] ** 2, 0.0))
    floor = max(sig0 / red, 1e-6)         # bruit attendu dans l'image lissée
    w = np.clip((std_r / floor - 2.0) / (5.0 - 2.0), 0.0, 1.0)    # bruit→0, signal→1
    w = gaussian_filter1d(w, 4.0)
    out[keep] = (rank * w[kr]).astype(np.float32)
    return out


def fnrgf(lum, cx, cy, lunar_r=0.0, order=8, smooth=8.0, denoise=2.0, rmax_px=0.0):
    """FNRGF — Fourier Normalizing Radial-Graded Filter (Druckmüllerová, Morgan &
    Habbal 2011). Moyenne ET écart-type locaux approchés par une série de Fourier
    azimutale d'ordre `order` → l'« attendu » suit les streamers : pas de lobes
    sombres. order faible = grandes structures (naturel) ; order élevé = très
    aplati. Retourne [0,1] ; disque + au-delà → sombre."""
    a = lum.astype(np.float32)
    H, W = a.shape
    yy, xx = np.indices((H, W)).astype(np.float32)
    r = np.hypot(xx - cx, yy - cy)
    th = np.arctan2(yy - cy, xx - cx)
    rmax = float(rmax_px) if rmax_px and rmax_px > 0 else float(r.max())
    nb = int(r.max()) + 1
    ri = np.clip(r.astype(np.int32), 0, nb - 1)
    keep = (r >= lunar_r) & (r <= rmax) & np.isfinite(a)
    order = max(1, int(order))

    idx = ri[keep]
    val = a[keep].astype(np.float64)
    ang = th[keep]
    cnt = np.bincount(idx, minlength=nb).astype(np.float64)
    good = cnt > 0
    allidx = np.arange(nb)

    def binmean(w):
        s = np.bincount(idx, weights=w, minlength=nb)
        m = np.zeros(nb)
        m[good] = s[good] / cnt[good]
        if good.sum() >= 2:
            m = np.interp(allidx, allidx[good], m[good])
        return gaussian_filter1d(m, max(1.0, smooth))

    # coefficients de Fourier (en azimut) par anneau, pour I et pour I²
    a0  = binmean(val)
    a02 = binmean(val * val)
    cs, sn, cs2, sn2 = [], [], [], []
    for k in range(1, order + 1):
        ck = np.cos(k * ang); sk = np.sin(k * ang)
        cs.append(2.0 * binmean(val * ck));  sn.append(2.0 * binmean(val * sk))
        cs2.append(2.0 * binmean(val * val * ck)); sn2.append(2.0 * binmean(val * val * sk))

    # reconstruction des séries au niveau de chaque pixel
    meanI  = a0[ri].copy()
    meanI2 = a02[ri].copy()
    for k in range(1, order + 1):
        ck = np.cos(k * th); sk = np.sin(k * th)
        meanI  += cs[k - 1][ri]  * ck + sn[k - 1][ri]  * sk
        meanI2 += cs2[k - 1][ri] * ck + sn2[k - 1][ri] * sk

    std = np.sqrt(np.maximum(meanI2 - meanI ** 2, 0.0))
    if denoise > 0:
        std = np.maximum(std, denoise * _noise_sigma(a))
    std = np.maximum(std, 1e-6)
    field = (a - meanI) / std
    v = field[keep]
    lo = float(np.percentile(v, 0.5)); hiq = float(np.percentile(v, 99.5))
    if hiq <= lo:
        hiq = lo + 1e-6
    out = np.clip((field - lo) / (hiq - lo), 0.0, 1.0)
    out[~keep] = 0.0
    return out.astype(np.float32)


# ── Détail tangentiel — ACHF simplifié (Druckmüller 2006) ─────────────────────────
def tangential_detail(lum, cx, cy, lunar_r=0.0, arc_deg=8.0, strength=1.0,
                      denoise=2.0, rmax_px=0.0):
    """Masque flou LE LONG des arcs de cercle, soustrait de l'image.

    En espace polaire (r, θ), on floute uniquement selon θ (arc de `arc_deg`
    degrés) : le flou suit les cercles, jamais les rayons. La différence
    I − flou_θ(I) ne contient donc QUE les structures radiales (streamers,
    jets, plumes) — le gradient radial est intact. Coring anti-bruit :
    |détail| < denoise·σ_n → 0. Seul le DÉTAIL est rééchantillonné
    (polaire aller-retour) ; l'image de base reste nette au pixel près.
    Retour : image + force × détail (échelle d'origine, ≥ 0)."""
    a = lum.astype(np.float32)
    H, W = a.shape
    yy, xx = np.indices((H, W)).astype(np.float32)
    r = np.hypot(xx - cx, yy - cy)
    maxR = float(r.max())
    nr = int(maxR) + 2
    na = int(min(4096, max(720, 2.0 * np.pi * maxR)))
    th = (2.0 * np.pi * np.arange(na, dtype=np.float32) / na)[:, None]
    rr = np.arange(nr, dtype=np.float32)[None, :]
    Y = cy + rr * np.sin(th)
    X = cx + rr * np.cos(th)
    P = map_coordinates(a, [Y, X], order=1, mode="nearest").astype(np.float32)
    sig_th = max(1.0, na * float(arc_deg) / 360.0)
    Pb = gaussian_filter1d(P, sig_th, axis=0, mode="wrap")   # flou angulaire pur
    d = P - Pb
    if denoise > 0:   # coring : sous le plancher de bruit → 0 (pas de grain remonté)
        floor = float(denoise) * _noise_sigma(a)
        d = np.sign(d) * np.maximum(np.abs(d) - floor, 0.0)
    # fenêtre radiale : rien sous le limbe, extinction douce au rayon max
    ramp = np.clip((rr - max(float(lunar_r), 0.0)) / 20.0, 0.0, 1.0)
    if rmax_px and rmax_px > 0:
        ramp = ramp * (1.0 - np.clip((rr - float(rmax_px)) / 20.0, 0.0, 1.0))
    d *= ramp
    # retour en cartésien : couture angulaire 0 = 2π par duplication de ligne
    dp = np.vstack([d, d[:1]])
    ang = np.arctan2(yy - cy, xx - cx)
    ang = np.where(ang < 0, ang + 2.0 * np.pi, ang) * (na / (2.0 * np.pi))
    det = map_coordinates(dp, [ang, r], order=1, mode="nearest").astype(np.float32)
    return np.maximum(a + float(strength) * det, 0.0)


# ── Image loader ──────────────────────────────────────────────────────────────
class ImageLoadWorker(QThread):
    done  = pyqtSignal(object, int, int)
    error = pyqtSignal(str)

    def __init__(self, siril, path=None):
        super().__init__()
        self.siril = siril
        self.path = path
        self.bitpix = None            # profondeur du FITS reçu (postérisation ?)

    def run(self):
        tmp_path = None
        try:
            if self.path:
                src = self.path
            else:
                tmp = tempfile.mktemp(prefix="corona_")
                tmp_path = src = tmp + ".fit"
                with self.siril.image_lock():
                    self.siril.cmd(f'save "{tmp.replace(os.sep, "/")}"')
            with fits.open(src) as hdul:
                self.bitpix = int(hdul[0].header.get("BITPIX", -32))
                raw = hdul[0].data.astype(np.float32)
            if raw.ndim == 3 and raw.shape[2] in (3, 4) and raw.shape[0] not in (3, 4):
                raw = raw.transpose(2, 0, 1)
            if raw.ndim == 3 and raw.shape[0] == 4:
                raw = raw[:3]
            mx = float(np.max(raw))
            if mx > 1.0:
                raw = raw / (65535.0 if mx <= 65535.0 else mx)
            if raw.ndim == 3 and raw.shape[0] == 3:
                if np.array_equal(raw[0], raw[1]) and np.array_equal(raw[1], raw[2]):
                    raw = raw[0]
            H, W = raw.shape[-2], raw.shape[-1]
            self.done.emit(raw, W, H)
        except Exception as e:
            self.error.emit(str(e))
        finally:
            if tmp_path and os.path.exists(tmp_path):
                try: os.remove(tmp_path)
                except Exception: pass


# ── MGN worker ───────────────────────────────────────────────────────────────
class MGNWorker(QThread):
    progress = pyqtSignal(int, str)
    done     = pyqtSignal(object)
    error    = pyqtSignal(str)

    def __init__(self, data, params, parent=None):
        super().__init__(parent)
        self.data = data
        self.p = params

    def run(self):
        try:
            p = self.p
            chans = to_channels(self.data)
            lum = chans[0] if len(chans) == 1 else np.mean(chans, axis=0).astype(np.float32)

            if p["method"] == "rhef":
                self.progress.emit(20, "RHEF — égalisation par anneau…")
                L = rhef(lum, p["cx"], p["cy"], p["lunar_r"], p["denoise"], p["rmax"])
                self.progress.emit(85, "RHEF — fini.")
            elif p["method"] == "fnrgf":
                self.progress.emit(20, "FNRGF — séries de Fourier azimutales…")
                L = fnrgf(lum, p["cx"], p["cy"], p["lunar_r"], p["order"],
                          p["smooth"], p["denoise"], p["rmax"])
                self.progress.emit(85, "FNRGF — fini.")
            elif p["method"] == "tangential":
                self.progress.emit(20, "Détail tangentiel — espace polaire…")
                L = tangential_detail(lum, p["cx"], p["cy"], p["lunar_r"],
                                      p["arc"], p["strength"], p["denoise"],
                                      p["rmax"])
                self.progress.emit(85, "Détail tangentiel — fini.")
            else:
                self.progress.emit(15, "MGN…")
                L = mgn(lum, p["scales"], p["k"], p["gamma"], p["h"], p["denoise"],
                        progress=lambda v, m: self.progress.emit(v, m))

            if len(chans) == 1 or not p["color"]:
                stacked = L
            else:
                # Couleur préservée : on applique la nouvelle luminance comme
                # facteur (la teinte R:G:B d'origine est conservée).
                ratio = np.clip(L / np.maximum(lum, 1e-6), 0.0, 20.0)
                stacked = np.stack([c * ratio for c in chans]).astype(np.float32)

            if (stacked.ndim == 3 and (p.get("neutral", 0.0) > 1e-6
                                     or p.get("neutral_wb", False))):
                self.progress.emit(90, "Neutralisation de la couronne…")
                stacked, g = neutralize_corona(stacked, p["neutral"],
                                               p.get("neutral_wb", True))
                self.progress.emit(91, f"Couronne neutralisée "
                                       f"(balance R/V/B ×{g[0]:.3f}/{g[1]:.3f}/{g[2]:.3f}).")

            self.progress.emit(92, "Normalisation…")
            if p["method"] == "tangential":
                # Tonalité du HDR PRÉSERVÉE : division par le max, aucune
                # remontée du point noir → le fond reste sombre et l'image
                # part telle quelle au GHS (seul le détail a été ajouté).
                hi = float(stacked.max())
                out = np.clip(stacked / max(hi, 1e-6), 0.0, 1.0).astype(np.float32)
            else:
                lo = float(np.percentile(stacked, 0.25))
                hi = float(np.percentile(stacked, 99.75))
                if hi <= lo:
                    hi = lo + 1e-6
                v = (stacked - lo) / (hi - lo)
                # La sortie de MGN/RHEF/FNRGF est déjà bornée : le cramage vient
                # UNIQUEMENT de cet étirement sur P99,75, qui écrête tout ce qui
                # dépasse — et ce qui dépasse est justement la couronne interne.
                # Coude doux (tanh) au lieu de l'écrêtage : le modelé est conservé.
                k = float(p.get("knee", 0.85))
                if k >= 0.999:
                    out = np.clip(v, 0.0, 1.0).astype(np.float32)
                else:
                    out = np.clip(np.where(v <= k, v,
                                           k + (1.0 - k) * np.tanh((v - k) / (1.0 - k))),
                                  0.0, 1.0).astype(np.float32)
            # ── Disque lunaire : noir pur ou gris du fond de ciel ─────────────
            # Les méthodes radiales forcent le disque à 0. Un noir absolu au
            # milieu d'un fond gris devient le point le plus contrasté de
            # l'image : l'œil y va au lieu d'aller à la couronne. Le remplir au
            # niveau du fond de ciel rend la Lune simplement silhouettée.
            lr = float(p.get("lunar_r", 0.0))
            if p.get("disc_grey", True) and lr > 4:
                H0, W0 = out.shape[-2], out.shape[-1]
                gy, gx = np.indices((H0, W0)).astype(np.float32)
                rr = np.hypot(gx - p["cx"], gy - p["cy"])
                lum_o = out if out.ndim == 2 else out.mean(axis=0)
                ref = (rr > lr * 1.05) & (lum_o > 0)
                if ref.any():
                    lvl = float(np.percentile(lum_o[ref], 20.0))   # fond de ciel
                    w = (1.0 - np.clip((rr - (lr - 8.0)) / 7.0, 0.0, 1.0)).astype(np.float32)
                    ww = w if out.ndim == 2 else w[None, :, :]
                    out = (out * (1.0 - ww) + lvl * ww).astype(np.float32)
                    self.progress.emit(97, f"Disque lunaire ramené au niveau du "
                                           f"fond de ciel ({lvl:.4f}).")
            self.progress.emit(100, "Terminé.")
            self.done.emit(out)
        except Exception as e:
            self.error.emit(str(e))


# ── Vue : clic-glissé gauche = centre (RHEF/FNRGF), molette = zoom, Maj+molette =
#         rayon lunaire, clic milieu glissé = déplacer la vue ──────────────────
class CenterView(QGraphicsView):
    center_at        = pyqtSignal(float, float)
    radius_delta     = pyqtSignal(int)
    radius_max_delta = pyqtSignal(int)

    def __init__(self, scene, parent=None):
        super().__init__(scene, parent)
        self.setTransformationAnchor(QGraphicsView.ViewportAnchor.AnchorUnderMouse)
        self.setResizeAnchor(QGraphicsView.ViewportAnchor.AnchorUnderMouse)
        self.setDragMode(QGraphicsView.DragMode.NoDrag)
        self.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform)
        self._pan = None
        self._dragging = False
        self.pick_enabled = False          # n'agit que pour les méthodes radiales

    def wheelEvent(self, event):
        up = event.angleDelta().y() > 0
        mods = event.modifiers()
        if self.pick_enabled and (mods & Qt.KeyboardModifier.ShiftModifier):
            self.radius_delta.emit(1 if up else -1)
        elif self.pick_enabled and (mods & Qt.KeyboardModifier.ControlModifier):
            self.radius_max_delta.emit(1 if up else -1)
        else:
            f = 1.15 if up else 1 / 1.15
            self.scale(f, f)

    def mousePressEvent(self, event):
        if self.pick_enabled and event.button() == Qt.MouseButton.LeftButton:
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


# ── Presets d'échelles ───────────────────────────────────────────────────────
SCALE_PRESETS = {
    "Fin (détail serré)":         (1.25, 2.5, 5, 10),
    "Standard":                   (1.25, 2.5, 5, 10, 20, 40),
    "Large (grandes structures)": (2.5, 5, 10, 20, 40, 80),
}


# ── Main window ────────────────────────────────────────────────────────────────
class CoronaWindow(QMainWindow):

    def __init__(self, siril):
        super().__init__()
        self.siril = siril
        self.setWindowTitle(f"Corona v1.1 (détail tangentiel) — v{VERSION}")
        self.setMinimumSize(980, 660)
        self.setStyleSheet(DARK_SS)

        self.data = None
        self.result = None
        self._loader = None
        self._worker = None
        self._pixitem = None
        self.cx = self.cy = 0.0       # centre (coords tableau), pour RHEF/FNRGF
        self._H = self._W = 0
        self._overlay = []

        self._build_ui()
        self._load_image()

    # =========================================================================
    def _build_ui(self):
        central = QWidget(); self.setCentralWidget(central)
        h = QHBoxLayout(central); h.setContentsMargins(6, 6, 6, 6); h.setSpacing(8)

        left = QWidget(); left.setFixedWidth(250)
        lv = QVBoxLayout(left); lv.setContentsMargins(0, 0, 0, 0); lv.setSpacing(6)

        row0 = QHBoxLayout()
        b_reload = QPushButton("⟳ Image Siril"); b_reload.clicked.connect(self._load_image)
        b_open = QPushButton("Fichier…"); b_open.clicked.connect(self._open_file)
        row0.addWidget(b_reload); row0.addWidget(b_open)
        lv.addLayout(row0)
        self.lbl_img = QLabel("—")
        self.lbl_img.setStyleSheet("font-size:9pt; color:#aaa;")
        lv.addWidget(self.lbl_img)

        mrow = QHBoxLayout()
        mrow.addWidget(QLabel("Méthode :"))
        self.cmb_method = QComboBox()
        self.cmb_method.addItems(["RHEF (auto, pour GHS)",
                                  "FNRGF (réglable)",
                                  "MGN (multi-échelle)",
                                  "Détail tangentiel (ACHF)"])
        self.cmb_method.currentIndexChanged.connect(self._on_method_changed)
        mrow.addWidget(self.cmb_method)
        lv.addLayout(mrow)

        # ── Groupe centre (méthodes radiales : RHEF / FNRGF) ──────────────────
        self.grp_center = QGroupBox("Centre — méthodes radiales (RHEF / FNRGF)")
        gn = QGridLayout(self.grp_center)
        b_auto = QPushButton("◎ Détecter le centre")
        b_auto.clicked.connect(self._auto_center)
        gn.addWidget(b_auto, 0, 0, 1, 2)
        self.lbl_center = QLabel("centre : —")
        self.lbl_center.setStyleSheet("font-size:9pt; color:#aaa;")
        gn.addWidget(self.lbl_center, 1, 0, 1, 2)
        gn.addWidget(QLabel("Rayon lunaire (px) :"), 2, 0)
        self.spn_lunar = QSpinBox(); self.spn_lunar.setRange(0, 20000)
        self.spn_lunar.valueChanged.connect(self._draw_overlay)
        gn.addWidget(self.spn_lunar, 2, 1)
        gn.addWidget(QLabel("Rayon max couronne (px) :"), 3, 0)
        self.spn_rmax = QSpinBox(); self.spn_rmax.setRange(0, 40000)
        self.spn_rmax.setSingleStep(10)
        self.spn_rmax.setToolTip("0 = pas de limite. Sinon, tout au-delà de ce rayon "
                                 "passe en noir (supprime les coins bruités).")
        self.spn_rmax.valueChanged.connect(self._draw_overlay)
        gn.addWidget(self.spn_rmax, 3, 1)
        lv.addWidget(self.grp_center)

        # ── Groupe FNRGF ──────────────────────────────────────────────────────
        self.grp_fnrgf = QGroupBox("FNRGF — séries de Fourier azimutales")
        gf = QGridLayout(self.grp_fnrgf)
        gf.addWidget(QLabel("Ordre Fourier :"), 0, 0)
        self.spn_order = QSpinBox(); self.spn_order.setRange(1, 40)
        self.spn_order.setValue(8)
        self.spn_order.setToolTip("Faible (2–4) = grandes structures, rendu naturel ; "
                                  "élevé (12–20) = très aplati, détail fin, supprime les lobes.")
        gf.addWidget(self.spn_order, 0, 1)
        gf.addWidget(QLabel("Lissage radial :"), 1, 0)
        self.spn_rsmooth = QDoubleSpinBox(); self.spn_rsmooth.setRange(1.0, 200.0)
        self.spn_rsmooth.setValue(8.0)
        gf.addWidget(self.spn_rsmooth, 1, 1)
        self.grp_fnrgf.setVisible(False)
        lv.addWidget(self.grp_fnrgf)

        # ── Groupe Détail tangentiel ──────────────────────────────────────────
        self.grp_tang = QGroupBox("Détail tangentiel — ACHF simplifié")
        gt = QGridLayout(self.grp_tang)
        gt.addWidget(QLabel("Arc de flou (°) :"), 0, 0)
        self.spn_arc = QDoubleSpinBox(); self.spn_arc.setRange(1.0, 45.0)
        self.spn_arc.setSingleStep(1.0); self.spn_arc.setValue(8.0)
        self.spn_arc.setToolTip("Longueur d'arc du flou circulaire. Petit = détail "
                                "fin ; grand = structures radiales larges.")
        gt.addWidget(self.spn_arc, 0, 1)
        gt.addWidget(QLabel("Force :"), 1, 0)
        self.spn_tstrength = QDoubleSpinBox(); self.spn_tstrength.setRange(0.1, 5.0)
        self.spn_tstrength.setSingleStep(0.1); self.spn_tstrength.setValue(1.0)
        gt.addWidget(self.spn_tstrength, 1, 1)
        self.grp_tang.setVisible(False)
        lv.addWidget(self.grp_tang)

        # ── Groupe MGN ────────────────────────────────────────────────────────
        self.grp_mgn = QGroupBox("MGN — Multi-Scale Gaussian Normalization")
        g = QGridLayout(self.grp_mgn)
        g.addWidget(QLabel("Échelles :"), 0, 0)
        self.cmb_scales = QComboBox()
        self.cmb_scales.addItems(list(SCALE_PRESETS.keys()))
        self.cmb_scales.setCurrentIndex(1)
        g.addWidget(self.cmb_scales, 0, 1)
        g.addWidget(QLabel("Contraste k :"), 1, 0)
        self.spn_k = QDoubleSpinBox(); self.spn_k.setRange(0.1, 3.0)
        self.spn_k.setSingleStep(0.1); self.spn_k.setValue(0.7)
        g.addWidget(self.spn_k, 1, 1)
        g.addWidget(QLabel("Gamma global :"), 2, 0)
        self.spn_gamma = QDoubleSpinBox(); self.spn_gamma.setRange(1.0, 6.0)
        self.spn_gamma.setSingleStep(0.1); self.spn_gamma.setValue(3.2)
        g.addWidget(self.spn_gamma, 2, 1)
        g.addWidget(QLabel("Global ↔ détail (h) :"), 3, 0)
        self.spn_h = QDoubleSpinBox(); self.spn_h.setRange(0.0, 1.0)
        self.spn_h.setSingleStep(0.05); self.spn_h.setValue(0.7)
        g.addWidget(self.spn_h, 3, 1)
        self.grp_mgn.setVisible(False)
        lv.addWidget(self.grp_mgn)

        # ── Commun ────────────────────────────────────────────────────────────
        grp_com = QGroupBox("Commun")
        gco = QGridLayout(grp_com)
        gco.addWidget(QLabel("Réduction du bruit :"), 0, 0)
        self.spn_denoise = QDoubleSpinBox(); self.spn_denoise.setRange(0.0, 10.0)
        self.spn_denoise.setSingleStep(0.5); self.spn_denoise.setValue(2.0)
        self.spn_denoise.setToolTip("Plancher de bruit (× σ estimé). 0 = aucun ; "
                                    "plus haut = moins de grain dans les zones plates.")
        gco.addWidget(self.spn_denoise, 0, 1)
        gco.addWidget(QLabel("Coude hautes lumières :"), 1, 0)
        self.spn_knee = QDoubleSpinBox(); self.spn_knee.setRange(0.50, 1.0)
        self.spn_knee.setSingleStep(0.05); self.spn_knee.setDecimals(2)
        self.spn_knee.setValue(0.85)
        self.spn_knee.setToolTip(
            "Compression douce des hautes lumières au lieu de l'écrêtage.\n"
            "Plus bas = couronne interne moins cramée (0,75 supprime\n"
            "totalement l'écrêtage). 1,00 = ancien comportement.\n"
            "Sans effet en mode Détail tangentiel (tonalité déjà préservée).")
        gco.addWidget(self.spn_knee, 1, 1)
        gco.addWidget(QLabel("Neutraliser la couronne :"), 2, 0)
        self.spn_neutral = QDoubleSpinBox(); self.spn_neutral.setRange(0.0, 1.0)
        self.spn_neutral.setSingleStep(0.1); self.spn_neutral.setValue(0.8)
        self.spn_neutral.setToolTip(
            "La couronne K est physiquement BLANCHE : sa teinte est un artefact,\n"
            "que Corona rend visible en remontant la couronne faible.\n"
            "La désaturer maximise le RAPPORT de saturation avec les\n"
            "protubérances rouges (Hα préservé) — la saturation absolue se\n"
            "remonte ensuite d'un curseur. 0 = inchangé · 1 = neutre.")
        gco.addWidget(self.spn_neutral, 2, 1)
        self.chk_wb = QCheckBox("Retirer la dominante (balance des blancs)")
        self.chk_wb.setChecked(True)
        self.chk_wb.setToolTip("Égalise les canaux sur la couronne, en excluant les "
                               "pixels rouges pour ne pas blanchir les protubérances.")
        gco.addWidget(self.chk_wb, 3, 0, 1, 2)
        self.chk_disc = QCheckBox("Disque lunaire au niveau du fond (pas noir)")
        self.chk_disc.setChecked(True)
        self.chk_disc.setToolTip(
            "Les méthodes radiales forcent le disque à noir pur. Au milieu d'un\n"
            "fond gris, ce noir absolu devient le point le plus contrasté de\n"
            "l'image et capte le regard au détriment de la couronne.\n"
            "Coché : le disque est ramené au niveau du fond de ciel — la Lune est\n"
            "simplement silhouettée, comme à l'œil nu.\n"
            "Décoché : noir pur (comportement d'origine).")
        gco.addWidget(self.chk_disc, 4, 0, 1, 2)
        self.chk_color = QCheckBox("Préserver la couleur")
        self.chk_color.setChecked(True)
        gco.addWidget(self.chk_color, 5, 0, 1, 2)
        lv.addWidget(grp_com)

        hint = QLabel("Centre = clic-glissé · Maj+molette = rayon lunaire · Ctrl+molette = "
                      "rayon max. RHEF = aplatir (sans réglage) → GHS. Tangentiel = affûter "
                      "les structures radiales (combo : RHEF, recharger, puis tangentiel).")
        hint.setWordWrap(True); hint.setStyleSheet("font-size:8pt; color:#888;")
        lv.addWidget(hint)

        lv.addStretch()
        self.btn_apply = QPushButton("▶  Appliquer"); self.btn_apply.setObjectName("BtnGo")
        self.btn_apply.clicked.connect(self._apply)
        lv.addWidget(self.btn_apply)
        self.btn_save = QPushButton("Enregistrer + charger dans Siril")
        self.btn_save.setEnabled(False); self.btn_save.clicked.connect(self._save_load)
        lv.addWidget(self.btn_save)
        h.addWidget(left)

        right = QWidget()
        rv = QVBoxLayout(right); rv.setContentsMargins(0, 0, 0, 0); rv.setSpacing(4)
        self._scene = QGraphicsScene()
        self._view = CenterView(self._scene)
        self._view.pick_enabled = True            # RHEF (radial) par défaut
        self._view.center_at.connect(self._on_center_at)
        self._view.radius_delta.connect(self._on_radius_delta)
        self._view.radius_max_delta.connect(self._on_radius_max_delta)
        rv.addWidget(self._view, stretch=1)
        self.prog = QProgressBar(); rv.addWidget(self.prog)
        self.log = QTextEdit(); self.log.setReadOnly(True); self.log.setFixedHeight(96)
        rv.addWidget(self.log)
        h.addWidget(right, stretch=1)

        self.statusBar().showMessage("Chargement…")

    # =========================================================================
    #  Image
    # =========================================================================
    def _load_image(self, path=None):
        self.statusBar().showMessage("Chargement de l'image…")
        self._loader = ImageLoadWorker(self.siril, path if isinstance(path, str) else None)
        self._loader.done.connect(self._on_image)
        self._loader.error.connect(lambda e: self._log(f"✗ {e}"))
        self._loader.start()

    def _open_file(self):
        f, _ = QFileDialog.getOpenFileName(self, "Ouvrir un FITS", os.getcwd(),
                                           "FITS (*.fit *.fits *.fts *.FIT *.FITS *.FTS)")
        if f:
            self._load_image(f)

    def _on_image(self, data, W, H):
        self.data = data
        self.result = None
        self._H, self._W = H, W
        self.btn_save.setEnabled(False)
        self._show(autostretch_display(data))
        kind = "couleur" if data.ndim == 3 else "mono"
        self.lbl_img.setText(f"{W}×{H} px · {kind}")
        self._log(f"Image {W}×{H} ({kind}) chargée.")
        bp = getattr(self._loader, "bitpix", None)
        if bp in (8, 16):
            self._log("⚠ Image reçue en entiers 16 bits : la couronne faible n'a "
                      "que peu de niveaux → postérisation possible. Passez Siril "
                      "en 32 bits (commande set32bits) puis rechargez, ou ouvrez "
                      "directement le FITS 32 bits de FusionHDR via « Fichier… ».")
        self._auto_center()

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
    #  Méthode & centre (NRGF)
    # =========================================================================
    def _on_method_changed(self, idx):
        # 0 = RHEF, 1 = FNRGF, 2 = MGN, 3 = Détail tangentiel
        is_radial = idx in (0, 1, 3)             # méthodes qui exigent le centre
        self.grp_center.setVisible(is_radial)
        self.grp_fnrgf.setVisible(idx == 1)
        self.grp_mgn.setVisible(idx == 2)
        self.grp_tang.setVisible(idx == 3)
        self._view.pick_enabled = is_radial
        self._draw_overlay()

    def _auto_center(self):
        if self.data is None:
            return
        lum = self.data if self.data.ndim == 2 else np.mean(self.data, axis=0)
        d = detect_disk(lum)
        if d is not None:
            self.cx, self.cy, r = d
            self.spn_lunar.blockSignals(True)
            self.spn_lunar.setValue(int(round(r)))
            self.spn_lunar.blockSignals(False)
            self._log(f"Centre détecté : ({self.cx:.0f}, {self.cy:.0f}), rayon {r:.0f} px.")
        else:
            self.cx, self.cy = self._W / 2.0, self._H / 2.0
            self._log("Disque non détecté — centre au milieu. Ajustez au clic-glissé.")
        self.lbl_center.setText(f"centre : ({self.cx:.0f}, {self.cy:.0f})")
        self._draw_overlay()

    def _on_center_at(self, sx, sy):
        if self.data is None:
            return
        self.cx = float(sx); self.cy = float(self._H - 1 - sy)
        self.lbl_center.setText(f"centre : ({self.cx:.0f}, {self.cy:.0f})")
        self._draw_overlay()

    def _on_radius_delta(self, direction):
        step = max(1, int(self.spn_lunar.value() * 0.03))
        self.spn_lunar.setValue(self.spn_lunar.value() + direction * step)

    def _on_radius_max_delta(self, direction):
        cur = self.spn_rmax.value()
        if cur == 0:                      # premier cran : démarre à une valeur visible
            if direction < 0:
                return
            cur = max(int(self.spn_lunar.value() * 2.5),
                      int(min(self._W, self._H) * 0.4), 50)
            self.spn_rmax.setValue(cur)
            return
        step = max(5, int(cur * 0.03))
        self.spn_rmax.setValue(max(0, cur + direction * step))

    def _draw_overlay(self):
        for it in self._overlay:
            self._scene.removeItem(it)
        self._overlay.clear()
        if self.data is None or self.cmb_method.currentIndex() not in (0, 1, 3):
            return
        sx = self.cx; sy = self._H - 1 - self.cy
        pen = QPen(QColor("#ff3333")); pen.setCosmetic(True); pen.setWidth(2)
        L = 12
        self._overlay.append(self._scene.addLine(sx - L, sy, sx + L, sy, pen))
        self._overlay.append(self._scene.addLine(sx, sy - L, sx, sy + L, pen))
        lr = self.spn_lunar.value()
        if lr > 0:
            penc = QPen(QColor("#33aaff")); penc.setCosmetic(True); penc.setWidth(2)
            self._overlay.append(self._scene.addEllipse(sx - lr, sy - lr, 2*lr, 2*lr, penc))
        rm = self.spn_rmax.value()
        if rm > 0:
            penm = QPen(QColor("#ffaa33")); penm.setCosmetic(True); penm.setWidth(2)
            penm.setStyle(Qt.PenStyle.DashLine)
            self._overlay.append(self._scene.addEllipse(sx - rm, sy - rm, 2*rm, 2*rm, penm))

    # =========================================================================
    #  Traitement
    # =========================================================================
    def _apply(self):
        if self.data is None:
            self._log("Aucune image.")
            return
        method = ("rhef", "fnrgf", "mgn", "tangential")[self.cmb_method.currentIndex()]
        params = {
            "method":  method,
            "denoise": self.spn_denoise.value(),
            "knee":    self.spn_knee.value(),
            "neutral":    self.spn_neutral.value(),
            "neutral_wb": self.chk_wb.isChecked(),
            "disc_grey":  self.chk_disc.isChecked(),
            "color":   self.chk_color.isChecked(),
            # MGN
            "scales": SCALE_PRESETS[self.cmb_scales.currentText()],
            "k":      self.spn_k.value(),
            "gamma":  self.spn_gamma.value(),
            "h":      self.spn_h.value(),
            # RHEF / FNRGF / tangentiel
            "cx": self.cx, "cy": self.cy,
            "lunar_r": float(self.spn_lunar.value()),
            "order":   self.spn_order.value(),
            "smooth":  self.spn_rsmooth.value(),
            "rmax":    float(self.spn_rmax.value()),
            # tangentiel
            "arc":      self.spn_arc.value(),
            "strength": self.spn_tstrength.value(),
        }
        self.btn_apply.setEnabled(False)
        self.btn_save.setEnabled(False)
        self.prog.setValue(0)
        if method == "rhef":
            self._log(f"━━ RHEF (centre {self.cx:.0f},{self.cy:.0f}, "
                      f"rayon {self.spn_lunar.value()})")
        elif method == "fnrgf":
            self._log(f"━━ FNRGF (centre {self.cx:.0f},{self.cy:.0f}, "
                      f"rayon {self.spn_lunar.value()}, ordre {self.spn_order.value()}, "
                      f"lissage {self.spn_rsmooth.value():.0f})")
        elif method == "tangential":
            self._log(f"━━ Détail tangentiel (centre {self.cx:.0f},{self.cy:.0f}, "
                      f"arc {self.spn_arc.value():.0f}°, "
                      f"force {self.spn_tstrength.value():.1f})")
        else:
            self._log(f"━━ MGN (échelles {self.cmb_scales.currentText()}, "
                      f"k={params['k']:.2g}, γ={params['gamma']:.2g}, h={params['h']:.2g})")
        self._worker = MGNWorker(self.data, params)
        self._worker.progress.connect(lambda v, m: (self.prog.setValue(v), self._log(m)))
        self._worker.done.connect(self._on_done)
        self._worker.error.connect(lambda e: (self._log(f"✗ {e}"),
                                              self.btn_apply.setEnabled(True)))
        self._worker.start()

    def _on_done(self, result):
        self.result = result
        self.btn_apply.setEnabled(True)
        self.btn_save.setEnabled(True)
        self._show(autostretch_display(result))
        self._log("✓ Traitement appliqué. Ajustez ou enregistrez (puis GHS dans Siril).")

    def _save_load(self):
        if self.result is None:
            return
        default = str(Path(os.getcwd()) / "corona.fit")   # répertoire de travail de Siril
        f, _ = QFileDialog.getSaveFileName(self, "Enregistrer", default, "FITS (*.fit *.fits)")
        if not f:
            return
        try:
            save_fits(self.result, f)
            self._log(f"✓ Enregistré : {f}")
            try:
                self.siril.cmd(f'load "{Path(f).with_suffix("")}"')   # guillemets : espaces
                self.siril.log(f"Corona — résultat chargé : {Path(f).name}")
                self._log("✓ Chargé dans Siril.")
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
    win = CoronaWindow(siril)
    win.show()
    app.exec()


if __name__ == "__main__":
    main()
