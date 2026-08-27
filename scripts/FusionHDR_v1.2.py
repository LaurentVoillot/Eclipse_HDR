##############################################
# FusionHDR v1.2 — zones non saturées (spatial)
# Fusion HDR radiométrique de poses bracketées
# Version 1.2.0 — variante expérimentale (basée sur v1.1)
##############################################
#
# SPDX-License-Identifier: GPL-3.0-or-later
#
# FusionHDR.py et FusionHDR_v1.1.py restent inchangés — cette variante
# ajoute le mode « Zones non saturées » (mis par défaut) :
#   1. masque DUR de saturation par pose (valeur ≥ seuil % pleine échelle) ;
#   2. masque ÉTENDU d'une marge (px) — la saturation contamine son
#      voisinage : bavure, halo de diffusion, haut de plage non linéaire ;
#   3. fondu SPATIAL (rampe sur la distance au masque, pas sur l'intensité)
#      → le relais entre poses se fait loin de la saturation, en pleine
#      zone linéaire où les poses concordent : pas de marche le long des
#      isophotes (anti-paliers structurel).
# Parmi les poses valides : pondération variance minimale (MLE, cf. v1.1),
# calibration croisée et coude doux conservés.
#
# Crédits
# -------
#   • Interface : PyQt6 / conventions VeraLux
#   • Méthode : fusion HDR radiométrique linéaire (échelle par temps de pose)
#   • Variance minimale : d'après Granados et al. 2010 — poids =
#     1/Var(radiance), modèle de bruit Var(I) = a + b·I estimé par pose

"""
FusionHDR — combine plusieurs poses de luminosités différentes (déjà
alignées et calibrées) en une seule image HDR linéaire 32 bits, chargée
dans Siril pour traitement ultérieur (courbes, etc.).

Usage général (pas seulement éclipse) :
  • Éclipse solaire — couronne + protubérances
  • Lunaire — clair de Terre, terminateur, éclipse de Lune
  • Planétaire — planète + satellites
  • Ciel profond — cœur brillant + extensions faibles (M42, M31…)

Principe
--------
Chaque pose est ramenée à une radiance commune (valeur / temps_de_pose),
puis combinée pixel à pixel :
  • Remplacement par seuil — la plus longue pose non saturée gagne
  • Mélange pondéré — moyenne pondérée par l'exposition (SNR), hors saturation
  • Variance minimale (défaut, v1.1) — poids = 1/Var(radiance) par pixel :
    l'estimateur statistiquement optimal. Le modèle de bruit Var(I)=a+b·I
    (lecture + photons) est mesuré automatiquement sur chaque pose ;
    aucun réglage. Passe continûment de w∝t² (zones faibles, bruit de
    lecture) à w∝t (zones brillantes, bruit de photons).

Les images doivent être **déjà recalées** (même cadrage) et linéaires.
L'alignement et la calibration restent à faire en amont (Siril, etc.).

Compatibilité
-------------
• Siril 1.3+   • Python 3.10+ (via sirilpy)
• Dépendances : numpy, astropy, PyQt6
"""

import sys
import os

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
from scipy.ndimage import (gaussian_filter, gaussian_filter1d, laplace, zoom,
                           distance_transform_edt, label, binary_closing)

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QGridLayout, QLabel, QPushButton, QComboBox, QCheckBox, QSpinBox,
    QDoubleSpinBox, QProgressBar, QTextEdit, QFileDialog, QGroupBox,
    QGraphicsView, QGraphicsScene, QGraphicsPixmapItem,
    QTableWidget, QTableWidgetItem, QHeaderView, QAbstractItemView
)
from PyQt6.QtCore import Qt, QThread, pyqtSignal, QSettings, QTimer
from PyQt6.QtGui import QImage, QPixmap, QPainter

VERSION = "1.2.0"

# ── Dark stylesheet ───────────────────────────────────────────────────────────
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
QComboBox          { background:#3c3c3c; border:1px solid #555; border-radius:3px;
                     padding:3px 6px; color:#d4d4d4; }
QComboBox QAbstractItemView { background:#2b2b2b; color:#d4d4d4;
                               selection-background-color:#555; }
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

# ── FITS helpers ────────────────────────────────────────────────────────────────
def load_fits(path: str) -> tuple[np.ndarray, fits.Header]:
    """Charge un FITS → float32, forme (H,W) ou (C,H,W). Valeurs d'origine."""
    with fits.open(path) as hdul:
        hdr  = hdul[0].header
        data = hdul[0].data.astype(np.float32)
    if data.ndim == 3 and data.shape[2] in (3, 4) and data.shape[0] not in (3, 4):
        data = data.transpose(2, 0, 1)          # (H,W,C) → (C,H,W)
    if data.ndim == 3 and data.shape[0] == 4:    # ignore alpha éventuel
        data = data[:3]
    return data, hdr

def get_exptime(hdr: fits.Header) -> Optional[float]:
    for kw in ("EXPTIME", "EXPOSURE", "EXP", "EXPOINUS"):
        if kw in hdr:
            try:
                v = float(hdr[kw])
                if v > 0:
                    return v
            except Exception:
                pass
    return None

def stack_info(hdr: fits.Header):
    """(n_empilées, livetime) si l'image est un empilement Siril, sinon (0, None).

    Siril ADDITIONNE les temps de pose lors d'un empilement (d'où l'option
    `-nosum`) : LIVETIME et parfois EXPTIME valent la somme. Or `stack med` ou
    `rej` produit une MOYENNE, pas une somme — le temps de pose à utiliser pour
    la fusion HDR est donc la pose UNITAIRE. Si les groupes de vitesses n'ont
    pas le même nombre d'images, une confusion ici décale les masters les uns
    par rapport aux autres et refabrique des bandes."""
    n = 0
    for kw in ("STACKCNT", "NCOMBINE"):
        if kw in hdr:
            try:
                n = int(float(hdr[kw]))
                break
            except Exception:
                pass
    lt = None
    if "LIVETIME" in hdr:
        try:
            lt = float(hdr["LIVETIME"])
        except Exception:
            pass
    return n, lt


def to_channels(data: np.ndarray) -> list[np.ndarray]:
    return [data] if data.ndim == 2 else [data[i] for i in range(data.shape[0])]

def max_chan(data: np.ndarray) -> np.ndarray:
    """Valeur max par pixel sur les canaux (pour le test de saturation)."""
    return data if data.ndim == 2 else data.max(axis=0)

def save_fits(data: np.ndarray, path: str, exptime=None, note=None):
    """Écrit le FITS 32 bits, avec EXPTIME si le temps de pose a un sens.

    EXPTIME n'est écrit que lorsque toutes les poses fusionnées partagent la
    même vitesse : le résultat est alors ramené à l'échelle ADU d'origine et
    se comporte comme une pose unique de cette durée — la fusion finale des
    masters lit donc directement la bonne valeur, sans saisie manuelle."""
    out = np.asarray(data, dtype=np.float32)
    hdu = fits.PrimaryHDU(out)
    if exptime is not None and exptime > 0:
        hdu.header["EXPTIME"] = (float(exptime), "s, ecrit par FusionHDR")
        hdu.header["EXPOSURE"] = (float(exptime), "s, alias de EXPTIME")
    if note:
        hdu.header["HISTORY"] = note[:70]
    hdu.writeto(path, overwrite=True)

# ── Affichage ────────────────────────────────────────────────────────────────
def _mtf(x, m, lo, hi):
    dist = hi - lo
    if dist < 1e-9:
        return np.where(x > lo, 1.0, 0.0).astype(np.float32)
    xp = np.clip((x - lo) / dist, 0.0, 1.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        y = ((m - 1.0) * xp) / ((2.0 * m - 1.0) * xp - m)
    return np.nan_to_num(y, nan=0.0, posinf=1.0, neginf=0.0).astype(np.float32)

def autostretch_display(data: np.ndarray) -> np.ndarray:
    """(H,W) ou (C,H,W), échelle quelconque → affichage [0,1] étiré."""
    chans = to_channels(data)
    mx = max(float(c.max()) for c in chans) or 1.0
    chans = [np.clip(c / mx, 0.0, 1.0) for c in chans]
    lum = chans[0] if len(chans) == 1 else np.mean(chans, axis=0)
    stride = max(1, lum.size // 500_000)
    samp = lum.flatten()[::stride]
    med = float(np.median(samp))
    mad = float(np.median(np.abs(samp - med))) * 1.4826 or 0.001
    c0  = max(0.0, med - 2.8 * mad)
    mt  = float(_mtf(np.float32(med - c0), 0.25, 0.0, 1.0))
    out = np.stack([_mtf(c, mt, c0, 1.0) for c in chans])
    return out[0] if len(chans) == 1 else out

def to_pixmap(disp: np.ndarray) -> QPixmap:
    """disp [0,1], (H,W) ou (C,H,W) en orientation FITS → QPixmap (Y vers le haut)."""
    if disp.ndim == 2:
        g = np.ascontiguousarray(np.flipud((np.clip(disp, 0, 1) * 255).astype(np.uint8)))
        h, w = g.shape
        img = QImage(g.tobytes(), w, h, w, QImage.Format.Format_Grayscale8)
    else:
        rgb = np.flipud((np.clip(disp, 0, 1) * 255).astype(np.uint8).transpose(1, 2, 0))
        rgb = np.ascontiguousarray(rgb[:, :, :3])
        h, w, _ = rgb.shape
        img = QImage(rgb.tobytes(), w, h, 3 * w, QImage.Format.Format_RGB888)
    return QPixmap.fromImage(img)

def smoothstep(lo: float, hi: float, x: np.ndarray) -> np.ndarray:
    if hi - lo < 1e-12:
        return (x >= hi).astype(np.float32)
    t = np.clip((x - lo) / (hi - lo), 0.0, 1.0)
    return (t * t * (3.0 - 2.0 * t)).astype(np.float32)


def parse_exposure(text):
    """Lit un temps de pose saisi à la main → secondes, ou None si invalide.

    Accepte la virgule décimale (clavier français), la notation fractionnaire
    naturelle « 1/250 », un « s » final et les espaces : 0.004 · 0,004 ·
    1/250 · 1/1000 s · 2s · 2 s.
    """
    if text is None:
        return None
    t = str(text).strip().lower().replace("s", "").replace(" ", "")
    t = t.replace(" ", "").replace(",", ".")
    if not t:
        return None
    try:
        if "/" in t:
            num, den = t.split("/", 1)
            num = float(num) if num else 1.0
            den = float(den)
            if den == 0:
                return None
            return num / den
        return float(t)
    except ValueError:
        return None


def _linfit(x, y):
    """Régression linéaire y ≈ a·x + b, CENTRÉE avant résolution.

    np.polyfit construit la matrice de Vandermonde sur les valeurs brutes :
    avec des radiances couvrant 4 décades (couronne interne vs externe), elle
    est mal conditionnée → « RankWarning » et pente potentiellement fausse.
    En centrant, le conditionnement devient trivial et le résultat exact.
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    if x.size < 2:
        return 1.0, 0.0
    xm = float(x.mean()); ym = float(y.mean())
    dx = x - xm
    den = float((dx * dx).sum())
    if not np.isfinite(den) or den <= 0.0:
        return 1.0, 0.0
    a = float((dx * (y - ym)).sum() / den)
    if not np.isfinite(a):
        return 1.0, 0.0
    return a, ym - a * xm


# ── Modèle de bruit (mode « Variance minimale ») ─────────────────────────────────
def estimate_noise_model(lum: np.ndarray, sat_level: float, tile: int = 64):
    """Modèle de bruit Var(I) ≈ a + b·I estimé de l'image elle-même.

    Courbe de transfert photonique « mono-image » : variance locale haute
    fréquence (MAD², robuste) vs intensité médiane, par tuiles ; on prend
    l'ENVELOPPE BASSE (quartile 25 %) par classe d'intensité pour ignorer la
    variance due aux vraies structures, puis on ajuste une droite.
    a ≈ bruit de lecture², b ≈ 1/gain (photons). Le facteur d'atténuation du
    filtre HF est identique pour toutes les poses → il se factorise et n'a
    aucun effet sur les poids relatifs.
    Retourne (a, b) dans l'échelle ADU² de l'image.
    """
    hp = lum - gaussian_filter(lum, 2.0)
    H, W = lum.shape
    t = max(16, int(tile))
    ms, vs = [], []
    for y0 in range(0, H - t + 1, t):
        for x0 in range(0, W - t + 1, t):
            tl = lum[y0:y0 + t, x0:x0 + t]
            m = float(np.median(tl))
            if m >= sat_level:            # tuile saturée : variance écrasée
                continue
            hh = hp[y0:y0 + t, x0:x0 + t]
            v = (1.4826 * float(np.median(np.abs(hh - np.median(hh))))) ** 2
            ms.append(m); vs.append(v)
    if len(ms) < 8:
        return 1e-12, 0.0
    ms = np.asarray(ms); vs = np.asarray(vs)
    # enveloppe basse par classe d'intensité (12 classes en quantiles)
    qs = np.quantile(ms, np.linspace(0.0, 1.0, 13))
    bx, by = [], []
    for i in range(12):
        sel = (ms >= qs[i]) & (ms <= qs[i + 1])
        if int(sel.sum()) >= 4:
            bx.append(float(np.median(ms[sel])))
            by.append(float(np.percentile(vs[sel], 25.0)))
    if len(bx) < 3:
        return max(float(np.percentile(vs, 25.0)), 1e-12), 0.0
    b, a = _linfit(bx, by)                                 # pente, ordonnée
    return max(float(a), 1e-12), max(float(b), 0.0)


# ── Détection du disque lunaire (mode Zones) — cf. Corona/EclipseAlign ───────────
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


def _prominence_detail(rad, rr, rri, nb, r_in, r_out, sigma_lum, k_sig=3.0):
    """Extrait les PROTUBÉRANCES d'une pose : excès localisés au-dessus de la
    ligne de base azimutale (moyenne par anneau, par canal) — le halo/anneau
    de couronne interne EST la ligne de base, donc il est retiré. Le résultat
    est fenêtré en douceur au-dessus de k_sig × σ (pas de bruit boosté).

    La bande est BORNÉE À LA HAUTEUR DES PROTUBÉRANCES au-dessus du limbe
    (r_in → r_out) et s'éteint progressivement à son bord extérieur : une
    bande trop large capterait les streamers de la couronne (qui sont eux
    aussi des excès azimutaux) et une coupure franche dessinerait un cercle.
    """
    chans = [rad] if rad.ndim == 2 else [rad[i] for i in range(rad.shape[0])]
    zone = (rr > r_in) & (rr < r_out)
    if not zone.any():
        return np.zeros_like(rad)
    idx = rri[zone]
    cnt = np.bincount(idx, minlength=nb).astype(np.float64)
    good = cnt > 0
    ai = np.arange(nb)
    dets = []
    for c in chans:
        s = np.bincount(idx, weights=c[zone].astype(np.float64), minlength=nb)
        m = np.zeros(nb)
        m[good] = s[good] / cnt[good]
        if good.sum() >= 2:
            m = np.interp(ai, ai[good], m[good])
        m = gaussian_filter1d(m, 3.0)
        d = np.zeros_like(c)
        d[zone] = np.maximum(c[zone] - m[idx], 0.0)
        dets.append(d)
    det = dets[0] if len(dets) == 1 else np.stack(dets)
    dl = det if det.ndim == 2 else det.mean(axis=0)
    t = np.clip((dl - k_sig * sigma_lum)
                / np.maximum(2.0 * sigma_lum, 1e-12), 0.0, 1.0)
    g = (t * t * (3.0 - 2.0 * t)).astype(np.float32)
    # Extinction douce sur les 40 % extérieurs de la bande : plus aucun bord
    # net possible, quel que soit le réglage de hauteur.
    r_fade = r_out - 0.4 * max(r_out - r_in, 1e-6)
    g = g * (1.0 - smoothstep(r_fade, r_out, rr))
    return det * (g if det.ndim == 2 else g[None, :, :])


def detect_moon_disk(lum):
    """Disque sombre enclos par la couronne → (cx, cy, r) ou None."""
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


# ── Fusion d'expositions de Mertens (2007) — pyramides laplaciennes ──────────────
#    Rendu tone-mappé direct (façon Photomatix « Naturel ») : pas de stade
#    linéaire → les hautes lumières (couronne interne) ne sont jamais cramées.
def _reduce(img):
    return gaussian_filter(img, 1.0, mode="nearest")[::2, ::2]

def _expand(img, shape):
    z = (shape[0] / img.shape[0], shape[1] / img.shape[1])
    return zoom(img, z, order=1, mode="nearest").astype(np.float32)

def _gauss_pyr(img, nlev):
    pyr = [img]
    for _ in range(nlev - 1):
        pyr.append(_reduce(pyr[-1]))
    return pyr

def _lap_pyr(img, nlev):
    g = _gauss_pyr(img, nlev)
    lap = [g[l] - _expand(g[l + 1], g[l].shape) for l in range(nlev - 1)]
    lap.append(g[-1])
    return lap

def _collapse(lap):
    out = lap[-1]
    for l in range(len(lap) - 2, -1, -1):
        out = lap[l] + _expand(out, lap[l].shape)
    return out

def _well_exposed(data, sigma):
    """Poids « bonne exposition » : max au gris moyen, → 0 vers noir ET blanc
    (c'est ce qui évite de cramer la couronne interne)."""
    e = np.exp(-((data - 0.5) ** 2) / (2.0 * sigma * sigma)).astype(np.float32)
    return e if data.ndim == 2 else np.prod(e, axis=0)


# ── HDR worker (fusion incrémentale, 1 image en mémoire à la fois) ──────────────
class HDRWorker(QThread):
    progress = pyqtSignal(int, str)
    done     = pyqtSignal(object, object)  # résultat, EXPTIME (ou None)
    error    = pyqtSignal(str)

    def __init__(self, params: dict, parent=None):
        super().__init__(parent)
        self.p = params
        self._abort = False

    def abort(self): self._abort = True

    def _radiance(self, data, exptime, bg_sub, bg_pct):
        if bg_sub:
            for c in to_channels(data):
                c -= np.percentile(c, bg_pct)
            np.clip(data, 0.0, None, out=data)
        return data / float(exptime)

    @staticmethod
    def _fit_pair(x, y):
        """Ajustement robuste y ≈ α·x + β (3 passes de sigma-clipping)."""
        sel = np.isfinite(x) & np.isfinite(y)
        alpha, beta = 1.0, 0.0
        for _ in range(3):
            if int(sel.sum()) < 100:
                break
            alpha, beta = _linfit(x[sel], y[sel])
            res = y - (alpha * x + beta)
            med = float(np.median(res[sel]))
            s = 1.4826 * float(np.median(np.abs(res[sel] - med)))
            if s <= 0:
                break
            sel = np.abs(res - med) < 3.0 * s
        return float(alpha), float(beta)

    def _cross_calibrate(self, items, fs, sat_lo, bg_sub, bg_pct):
        """Calibration croisée anti-bandes : corrections (gain, offset) par pose.

        Les zones de fusion changent de pose dominante le long des isophotes ;
        si les radiances I/t ne concordent pas exactement (EXIF arrondi,
        non-linéarité près de la saturation, fond de ciel variable), chaque
        transition laisse une MARCHE concentrique. Ici, chaque pose est
        régressée (robuste) sur la suivante plus longue, sur leurs pixels
        communs valides, puis les corrections sont chaînées vers la pose la
        plus longue (référence — meilleur SNR dans la couronne faible).
        Retourne {index: (A, B)} tel que radiance_corrigée = A·radiance + B."""
        n = len(items)
        order = sorted(range(n), key=lambda i: items[i][1])   # court → long
        corr = {i: (1.0, 0.0) for i in range(n)}
        lo_floor = 0.02 * fs                                  # au-dessus du bruit

        # ── Passe 1 : ajuster chaque paire voisine, SANS rien appliquer ──────
        fits = []                                             # (i_short, α, β, noms)
        for pos in range(len(order) - 1, 0, -1):
            i_long, i_short = order[pos], order[pos - 1]
            p_l, e_l = items[i_long]; p_s, e_s = items[i_short]
            d_l, _ = load_fits(p_l); d_s, _ = load_fits(p_s)
            v_l = max_chan(d_l); v_s = max_chan(d_s)
            valid = ((v_l > lo_floor) & (v_l < sat_lo) &
                     (v_s > lo_floor) & (v_s < sat_lo))
            r_l = self._radiance(d_l, e_l, bg_sub, bg_pct)
            r_s = self._radiance(d_s, e_s, bg_sub, bg_pct)
            rl = max_chan(r_l)[valid].ravel()
            rs = max_chan(r_s)[valid].ravel()
            if rl.size > 200_000:                             # sous-échantillonnage
                st = rl.size // 200_000
                rl = rl[::st]; rs = rs[::st]
            if rl.size < 500:
                self.progress.emit(2, f"⚠ Calibration : recouvrement insuffisant "
                                      f"{Path(p_s).name} ↔ {Path(p_l).name} — ignorée.")
                fits.append((i_short, 1.0, 0.0, Path(p_s).name, Path(p_l).name))
                continue
            alpha, beta = self._fit_pair(rl, rs)              # court ≈ α·long + β
            if not np.isfinite(alpha) or not np.isfinite(beta) or alpha <= 0:
                alpha, beta = 1.0, 0.0
            fits.append((i_short, alpha, beta, Path(p_s).name, Path(p_l).name))

        # ── Contrôle de cohérence global ─────────────────────────────────────
        # Cette calibration corrige des écarts RÉSIDUELS (quelques %) : EXIF
        # arrondis, non-linéarité, fond de ciel. Si les écarts atteignent un
        # facteur ~2 par cran, ce ne sont pas des résidus — les temps de pose
        # eux-mêmes sont faux. Les chaîner amplifierait la pose la plus courte
        # d'un facteur énorme (×177 sur 8 crans) et brûlerait toute l'image.
        # Dans ce cas : on n'applique RIEN (comportement homogène) et on le dit.
        gains = [1.0 / a for _, a, _, _, _ in fits if a > 0]
        gmed = float(np.median(gains)) if gains else 1.0
        if not (0.8 < gmed < 1.25):
            self.progress.emit(2,
                f"⚠ Calibration croisée DÉSACTIVÉE : écart médian de ×{gmed:.2f} par "
                f"cran entre poses voisines. Ce n'est pas un résidu — vos temps de "
                f"pose sont faux (rapport ~×{gmed:.1f} par cran). Saisissez les vraies "
                f"valeurs EXIF dans la colonne « Expo (s) » plutôt que de les estimer.")
            return corr

        # ── Passe 2 : chaînage vers la pose la plus longue (référence) ───────
        A_next, B_next = 1.0, 0.0
        for i_short, alpha, beta, nm_s, nm_l in fits:
            A = A_next / alpha
            B = B_next - A_next * beta / alpha
            corr[i_short] = (A, B)
            self.progress.emit(2, f"Calibration {nm_s} ↔ {nm_l} : "
                                  f"gain {1.0/alpha:.4f}, offset {-beta/alpha:.3g}")
            A_next, B_next = A, B
        return corr

    def run(self):
        try:
            items   = self.p["items"]          # [(path, exptime), …]
            mode    = self.p["mode"]   # "threshold" | "weighted" | "mertens" | "mle"
            fs      = float(self.p["fs"])
            sat_hi  = self.p["sat_frac"] * fs
            sat_lo  = max(0.0, (self.p["sat_frac"] - self.p["feather"]) * fs)
            bg_sub  = self.p["bg_sub"]
            bg_pct  = self.p["bg_pct"]
            norm    = self.p["normalize"]
            n = len(items)
            if n < 2:
                self.error.emit("Ajouter au moins deux poses.")
                return

            # Poses en entiers 8/16 bits ? La couronne faible (~1e-3 de la
            # pleine échelle) n'y a que ~65 niveaux → postérisation (terrasses
            # le long des isophotes) impossible à rattraper en aval.
            try:
                bps = {int(fits.getheader(p).get("BITPIX", -32)) for p, _ in items}
                if any(b in (8, 16) for b in bps):
                    self.progress.emit(1,
                        "⚠ Poses en entiers 8/16 bits détectées. Passez Siril en "
                        "32 bits (commande set32bits) et refaites conversion + "
                        "empilement des masters — sinon la postérisation résiduelle "
                        "persistera malgré la fusion 32 bits.")
            except Exception:
                pass

            # Calibration croisée anti-bandes (modes radiométriques seulement)
            corr = {i: (1.0, 0.0) for i in range(n)}
            if self.p.get("xcal", True) and mode != "mertens":
                self.progress.emit(2, "Calibration croisée des poses (anti-bandes)…")
                corr = self._cross_calibrate(items, fs, sat_lo, bg_sub, bg_pct)

            if mode == "threshold":
                # Tri par exposition croissante : on part de la plus courte
                # (la moins saturée) et chaque pose plus longue remplace là où
                # elle n'est pas saturée → la plus longue non saturée gagne.
                order = sorted(range(n), key=lambda i: items[i][1])
                p0, e0 = items[order[0]]
                self.progress.emit(int(100*1/(n+1)), f"Base : {Path(p0).name}")
                base, _ = load_fits(p0)
                hdr = self._radiance(base, e0, bg_sub, bg_pct)
                A0, B0 = corr[order[0]]
                hdr = A0 * hdr + B0
                for k, idx in enumerate(order[1:], start=2):
                    if self._abort:
                        return
                    path, exp = items[idx]
                    data, _ = load_fits(path)
                    if data.shape != hdr.shape:
                        self.error.emit(f"Dimensions différentes : {Path(path).name}")
                        return
                    val = max_chan(data)                       # pré-radiance
                    w   = 1.0 - smoothstep(sat_lo, sat_hi, val)  # 1 si OK, 0 si saturé
                    rad = self._radiance(data, exp, bg_sub, bg_pct)
                    Ai, Bi = corr[idx]
                    rad = Ai * rad + Bi
                    if hdr.ndim == 3:
                        w = w[None, :, :]
                    hdr = w * rad + (1.0 - w) * hdr
                    self.progress.emit(int(100*(k)/(n+1)), f"{k}/{n}  —  {Path(path).name}")

            elif mode == "mertens":
                # Fusion d'expositions (Mertens 2007) : pondération « bonne
                # exposition » × contraste, mélange multi-échelle (pyramides
                # laplaciennes). N'utilise NI les temps de pose NI le seuil de
                # saturation → les hautes lumières ne sont jamais cramées.
                sigma = 0.2
                paths = [it[0] for it in items]
                d0, _ = load_fits(paths[0]); d0 = d0 / fs
                H, W = d0.shape[-2], d0.shape[-1]
                nchan = 1 if d0.ndim == 2 else d0.shape[0]
                # Pyramide quasi complète : plus de niveaux = mélange plus complet
                # = moins de halo sombre au limbe (bord très contrasté).
                nlev = max(3, min(11, int(np.log2(min(H, W)))))
                del d0
                # Passe 1 : cartes de poids + somme
                wmaps = []; wsum = None
                for j, p in enumerate(paths):
                    if self._abort:
                        return
                    data, _ = load_fits(p)
                    if data.shape[-2:] != (H, W):
                        self.error.emit(f"Dimensions différentes : {Path(p).name}")
                        return
                    data = np.clip(data / fs, 0.0, 1.0)
                    lum = data if data.ndim == 2 else data.mean(axis=0)
                    contrast = np.abs(laplace(lum.astype(np.float32)))
                    wmap = ((contrast + 1e-4) * _well_exposed(data, sigma)
                            + 1e-12).astype(np.float32)
                    wmaps.append(wmap)
                    wsum = wmap if wsum is None else wsum + wmap
                    self.progress.emit(5 + int(35*(j+1)/n), f"poids {j+1}/{n}")
                wsum = np.where(wsum <= 0, 1e-12, wsum)
                # Passe 2 : pyramide laplacienne mélangée
                blend = [[None] * nchan for _ in range(nlev)]
                for j, p in enumerate(paths):
                    if self._abort:
                        return
                    data, _ = load_fits(p); data = np.clip(data / fs, 0.0, 1.0)
                    gw = _gauss_pyr(wmaps[j] / wsum, nlev)
                    chans = to_channels(data)
                    for ci in range(nchan):
                        lp = _lap_pyr(chans[ci].astype(np.float32), nlev)
                        for l in range(nlev):
                            contrib = gw[l] * lp[l]
                            blend[l][ci] = (contrib if blend[l][ci] is None
                                            else blend[l][ci] + contrib)
                    self.progress.emit(40 + int(55*(j+1)/n), f"fusion {j+1}/{n}")
                outs = [_collapse([blend[l][ci] for l in range(nlev)])
                        for ci in range(nchan)]
                hdr = np.clip(outs[0] if nchan == 1 else np.stack(outs), 0.0, 1.0)

            elif mode == "zones":
                # Zones non saturées : masque dur ≥ seuil, étendu d'une marge
                # (la saturation contamine son voisinage : bavure, halo,
                # non-linéarité du haut de plage), puis rampe de poids sur la
                # DISTANCE au masque. Le relais entre poses se fait ainsi loin
                # de la saturation, en zone linéaire → pas de marche le long
                # des isophotes. Pondération MLE parmi les poses valides.
                margin = float(self.p.get("zone_margin", 20))
                sfeather = max(float(self.p.get("zone_feather", 20)), 1e-6)
                short_i = min(range(n), key=lambda i: items[i][1])

                # ── Passe 1 : cercle lunaire de CHAQUE pose. La Lune peut
                # dériver entre poses (alignement non fait sur elle, ou
                # résidu) : ses pixels sombres — parfaitement « valides » et
                # sur-pondérés (val≈0 → 1/variance max) — écraseraient la
                # moyenne là où d'autres poses voient une protubérance →
                # taches noires au limbe. Règle : un pixel « Lune » d'une
                # pose ne vote pas ; là où TOUTES les poses voient la Lune
                # (intersection), on fusionne normalement (disque cohérent).
                circles = []
                for k, (path, exp) in enumerate(items, start=1):
                    if self._abort:
                        return
                    d0, _ = load_fits(path)
                    det = detect_moon_disk(max_chan(d0))
                    circles.append(det)
                    del d0
                    self.progress.emit(int(20*k/n),
                        f"Lune {k}/{n} : " + (f"r={det[2]:.0f}px" if det else
                        "non détectée (pose conservée, cercle médian appliqué)"))
                rs_det = [c[2] for c in circles if c]
                circles_prom = list(circles)      # cercles pour les PROTUBÉRANCES
                if rs_det:
                    r_med = float(np.median(rs_det))
                    circles = [c if (c and abs(c[2] - r_med) < 0.2 * r_med) else None
                               for c in circles]
                    # Cercle médian pour les poses non détectées — mais UNIQUEMENT
                    # pour l'extraction des protubérances, qui a juste besoin d'un
                    # centre et dont l'apport est additif (donc sans risque).
                    # PAS pour le masque lunaire : celui-ci arbitre le DÉSACCORD
                    # entre poses ; donner un cercle à une pose qui n'en revendique
                    # pas change qui vote à l'intérieur du disque, et donc le
                    # niveau du disque lui-même.
                    ok = [c for c in circles if c]
                    if ok and len(ok) < len(circles):
                        cx_m = float(np.median([c[0] for c in ok]))
                        cy_m = float(np.median([c[1] for c in ok]))
                        r_m = float(np.median([c[2] for c in ok]))
                        n_fill = sum(1 for c in circles if c is None)
                        circles_prom = [c if c else (cx_m, cy_m, r_m) for c in circles]
                        self.progress.emit(20,
                            f"Lune : cercle médian ({cx_m:.0f},{cy_m:.0f}) r={r_m:.0f}px "
                            f"appliqué à {n_fill} pose(s) pour les protubérances "
                            f"seulement (masque lunaire inchangé).")
                if not any(circles):
                    self.progress.emit(20, "⚠ Lune détectée sur aucune pose : masque "
                                           "lunaire désactivé.")
                H0 = W0 = None
                moon_grid = None
                m_list = None          # masque « je vois la Lune » par pose
                all_moon = None        # min des m_k : toutes les poses d'accord
                limbM = None
                limb_on = bool(self.p.get("limb_mono", False))
                prom_on = bool(self.p.get("prom_overlay", True))
                prom_n = max(1, min(int(self.p.get("prom_n", 3)), n))
                prom_h = max(0.02, float(self.p.get("prom_h", 0.15)))
                order_exp = sorted(range(n), key=lambda i: items[i][1])
                short_set = set(order_exp[:prom_n])   # les N poses les plus courtes
                detail_union = None            # union (max) des protubérances extraites

                num = den = None
                den2d = pot2d = rmax = None    # repli « tout saturé »
                has_allowed = None             # ≥1 pose non masquée à ce pixel
                for k, (path, exp) in enumerate(items, start=1):
                    if self._abort:
                        return
                    data, _ = load_fits(path)
                    if num is not None and data.shape != num.shape:
                        self.error.emit(f"Dimensions différentes : {Path(path).name}")
                        return
                    val = max_chan(data)                   # pré-radiance (ADU)
                    satm = val >= (self.p["sat_frac"] * fs)
                    if satm.any():
                        dist = distance_transform_edt(~satm)
                        wsp = np.clip((dist - margin) / sfeather,
                                      0.0, 1.0).astype(np.float32)
                    else:
                        wsp = np.ones(val.shape, np.float32)
                    # Masques lunaires : préparés une fois pour toutes les poses.
                    if moon_grid is None:
                        H0, W0 = val.shape
                        gy, gx = np.indices((H0, W0)).astype(np.float32)
                        moon_grid = (gy, gx)
                        m_list = []
                        for c in circles:
                            if c is None:            # cercle non détecté : cette
                                m_list.append(np.zeros((H0, W0), np.float32))
                                continue             # pose ne revendique rien
                            dd = np.hypot(gx - c[0], gy - c[1])
                            # 1 dans la Lune → 0 dehors, bord adouci sur 6 px, et
                            # décalé de +4 px : le rayon détecté (√(aire/π))
                            # sous-estime un peu le bord réel.
                            m_list.append(1.0 - np.clip((dd - (c[2] + 4.0)) / 6.0,
                                                        0.0, 1.0).astype(np.float32))
                        all_moon = m_list[0].copy()
                        for m in m_list[1:]:
                            all_moon = np.minimum(all_moon, m)
                        limbM = np.zeros((H0, W0), np.float32)
                        if limb_on and any(c is not None for c in circles):
                            mband = np.zeros((H0, W0), bool)
                            for c in circles:
                                if c is not None:
                                    dd = np.hypot(gx - c[0], gy - c[1])
                                    mband |= np.abs(dd - c[2]) < (margin + 8.0)
                            de = distance_transform_edt(~mband)
                            limbM = np.clip(1.0 - de / 10.0, 0.0, 1.0).astype(np.float32)
                            self.progress.emit(22, "Limbe mono-pose : bande de "
                                f"{int(mband.sum())} px prise de la pose la plus courte.")
                    # Masque lunaire CONTINU, fondé sur le DÉSACCORD entre poses :
                    #   • toutes les poses voient la Lune  → toutes votent
                    #     (disque cohérent et sombre, aucune zone sans poids) ;
                    #   • cette pose voit la Lune mais une autre voit la couronne
                    #     → cette pose se tait (c'est elle qui créait les taches) ;
                    #   • cette pose ne voit pas la Lune → elle vote pleinement.
                    m_k = m_list[k - 1]
                    allow = (1.0 - m_k + m_k * all_moon).astype(np.float32)
                    if limb_on and limbM is not None and (k - 1) != short_i:
                        allow = allow * (1.0 - limbM)
                    a_i, b_i = estimate_noise_model(val, sat_lo)
                    var = a_i + b_i * np.maximum(val, 0.0)
                    w2 = ((float(exp) ** 2) / np.maximum(var, 1e-12)).astype(np.float32)
                    w = w2 * wsp * allow
                    rad = self._radiance(data, exp, bg_sub, bg_pct)
                    Ai, Bi = corr[k - 1]
                    rad = Ai * rad + Bi
                    # Repli là où TOUT est saturé : le maximum des radiances
                    # calibrées. Pour un pixel saturé, valeur/t est une borne
                    # inférieure, et la plus serrée vient de la pose la plus
                    # courte — or au seuil de saturation cette borne ÉGALE la
                    # vraie valeur : le max se raccorde donc en continu au
                    # régime normal (une moyenne pondérée, elle, privilégierait
                    # la pose longue = borne lâche → marche visible).
                    # rmax ne retient que les poses autorisées ici : sinon, dans
                    # le disque lunaire, la lueur résiduelle divisée par une pose
                    # très courte remonte et dessine un anneau gris.
                    pw = w2 * allow
                    den2d = w.copy() if den2d is None else den2d + w
                    pot2d = pw.copy() if pot2d is None else pot2d + pw
                    am = allow > 0.01
                    amb = am if rad.ndim == 2 else am[None, :, :]
                    cand = np.where(amb, rad, -np.inf)
                    rmax = cand if rmax is None else np.maximum(rmax, cand)
                    has_allowed = am.copy() if has_allowed is None else (has_allowed | am)
                    if prom_on and (k - 1) in short_set:
                        c_k = circles_prom[k - 1]   # cercle comblé (voir passe 1)
                        if c_k is None:
                            self.progress.emit(int(100*k/(n+1)),
                                f"⚠ Protubérances : cercle lunaire non détecté sur "
                                f"{Path(path).name} — pose ignorée pour l'union.")
                        else:
                            rr = np.hypot(moon_grid[1] - c_k[0], moon_grid[0] - c_k[1])
                            nb_r = int(rr.max()) + 1
                            rri = np.clip(rr.astype(np.int32), 0, nb_r - 1)
                            r_in = c_k[2] - 4.0
                            r_out = c_k[2] * (1.0 + prom_h)   # hauteur physique
                            sig = (np.sqrt(np.maximum(var, 0.0))
                                   / max(float(exp), 1e-9)).astype(np.float32)
                            det = _prominence_detail(rad, rr, rri, nb_r,
                                                     r_in, r_out, sig, 3.0)
                            detail_union = det if detail_union is None \
                                else np.maximum(detail_union, det)
                    if (k - 1) == short_i:
                        w = np.maximum(w, 1e-12)   # plancher où tout sature
                    if data.ndim == 3:
                        w = w[None, :, :]
                    contrib_n = w * rad
                    if num is None:
                        num = contrib_n
                        den = np.broadcast_to(w, rad.shape).astype(np.float32).copy()
                    else:
                        num += contrib_n
                        den += w
                    pc_sat = 100.0 * float(satm.mean())
                    self.progress.emit(int(100*k/(n+1)),
                                       f"{k}/{n}  —  {Path(path).name}  "
                                       f"[zone saturée {pc_sat:.1f} % + marge {margin:.0f}px]")
                den = np.where(den <= 0, 1e-12, den)
                hdr = num / den
                # Comblement là où la saturation a tout retiré (et seulement
                # là où une pose voyait autre chose que la Lune, sinon anneau
                # gris autour du disque).
                ok2d = has_allowed if has_allowed is not None \
                    else np.ones(hdr.shape[-2:], bool)
                rmax = np.where(np.isfinite(rmax), rmax, 0.0)
                rho = den2d / np.maximum(pot2d, 1e-12)
                fb = np.clip(1.0 - rho / 0.05, 0.0, 1.0) * ok2d
                if fb.max() > 0:
                    n_px = int((fb > 0.5).sum())
                    if hdr.ndim == 3:
                        fb = fb[None, :, :]
                    hdr = (1.0 - fb) * hdr + fb * np.maximum(rmax, 0.0)
                    self.progress.emit(int(100*n/(n+1)),
                        f"Zones saturées partout : {n_px} px sur la borne "
                        f"inférieure (pose la plus courte).")
                # ── Superposition des protubérances SANS halo (proposition
                # utilisateur) : de chaque pose courte on n'extrait QUE les
                # excès localisés au-dessus de sa ligne de base azimutale (=
                # ses protubérances ; le halo/anneau EST la ligne de base,
                # donc retiré), puis UNION (max) des N époques — toutes les
                # protubérances apparaissent, aucune n'est érodée par un vote,
                # et les halos d'époques différentes ne se mélangent jamais.
                if prom_on and detail_union is not None:
                    du_lum = detail_union if detail_union.ndim == 2 \
                        else detail_union.mean(axis=0)
                    n_boost = int((du_lum > 0).sum())
                    if n_boost:
                        hdr = hdr + detail_union
                        self.progress.emit(int(100*n/(n+1)),
                            f"Protubérances : union de {prom_n} pose(s) courte(s), "
                            f"{n_boost} px ajoutés (sans halo).")

            elif mode == "mle":
                # Variance minimale (Granados 2010) : poids pixel = 1/Var(radiance)
                # = t² / (a + b·I). Le modèle (a, b) est mesuré sur chaque pose.
                # Optimal continûment : w∝t² en zone faible (bruit de lecture),
                # w∝t en zone brillante (bruit de photons). Sans réglage.
                short_i = min(range(n), key=lambda i: items[i][1])
                num = den = None
                for k, (path, exp) in enumerate(items, start=1):
                    if self._abort:
                        return
                    data, _ = load_fits(path)
                    if num is not None and data.shape != num.shape:
                        self.error.emit(f"Dimensions différentes : {Path(path).name}")
                        return
                    val = max_chan(data)                       # pré-radiance (ADU)
                    a_i, b_i = estimate_noise_model(val, sat_lo)
                    var = a_i + b_i * np.maximum(val, 0.0)
                    # Poids en CLOCHE (triangle de Debevec) au lieu du seuil :
                    # chaque pixel est un mélange de plusieurs poses sur TOUTE
                    # la dynamique → les désaccords résiduels entre poses se
                    # fondent en dégradés au lieu de faire des marches
                    # localisées (anti-paliers), et la pose la plus adaptée
                    # domine naturellement (bruit réduit en zone faible).
                    hat = (np.maximum(val, 0.0)
                           * np.maximum(1.0 - val / max(sat_hi, 1e-9), 0.0))
                    w = ((float(exp) ** 2) / np.maximum(var, 1e-12)) * hat
                    if (k - 1) == short_i:
                        w = np.maximum(w, 1e-12)   # anti-trou où tout sature
                    rad = self._radiance(data, exp, bg_sub, bg_pct)
                    Ai, Bi = corr[k - 1]
                    rad = Ai * rad + Bi
                    if data.ndim == 3:
                        w = w[None, :, :]
                    contrib_n = w * rad
                    if num is None:
                        num = contrib_n
                        den = np.broadcast_to(w, rad.shape).astype(np.float32).copy()
                    else:
                        num += contrib_n
                        den += w
                    self.progress.emit(int(100*k/(n+1)),
                                       f"{k}/{n}  —  {Path(path).name}  "
                                       f"[σ² = {a_i:.3g} + {b_i:.3g}·I]")
                den = np.where(den <= 0, 1e-12, den)
                hdr = num / den

            else:  # weighted — moyenne pondérée par l'exposition, hors saturation
                # La pose la plus courte reçoit un poids plancher : là où TOUTES
                # les poses saturent, le résultat retombe sur elle (pas de trou noir).
                short_i = min(range(n), key=lambda i: items[i][1])
                num = den = None
                for k, (path, exp) in enumerate(items, start=1):
                    if self._abort:
                        return
                    data, _ = load_fits(path)
                    if num is not None and data.shape != num.shape:
                        self.error.emit(f"Dimensions différentes : {Path(path).name}")
                        return
                    val = max_chan(data)
                    w   = float(exp) * (1.0 - smoothstep(sat_lo, sat_hi, val))
                    if (k - 1) == short_i:
                        w = np.maximum(w, 1e-6)
                    rad = self._radiance(data, exp, bg_sub, bg_pct)
                    Ai, Bi = corr[k - 1]
                    rad = Ai * rad + Bi
                    if data.ndim == 3:
                        w = w[None, :, :]
                    contrib_n = w * rad
                    if num is None:
                        num = contrib_n
                        den = np.broadcast_to(w, rad.shape).astype(np.float32).copy()
                    else:
                        num += contrib_n
                        den += w
                    self.progress.emit(int(100*k/(n+1)), f"{k}/{n}  —  {Path(path).name}")
                den = np.where(den <= 0, 1e-12, den)
                hdr = num / den

            # Normalisation de sortie — par un haut percentile (robuste aux
            # pixels chauds / pixels saturés jusque dans la pose la plus courte,
            # qui sinon, via le max brut, assombriraient toute l'image).
            # (Sauté pour Mertens, dont la sortie est déjà équilibrée dans [0,1].)
            #
            # v1.1 : au-dessus du COUDE, compression douce (tanh) au lieu de
            # l'écrêtage dur. La couronne interne dépasse largement 0,01 % des
            # pixels : l'ancien clip la coupait net (cœur cramé) alors que le
            # modelé existait. knee=1.0 → écrêtage (ancien comportement).
            # ── Poses toutes à la MÊME vitesse ────────────────────────────────
            # Le résultat n'est alors pas un « HDR » mais une simple combinaison.
            # On le remet à l'échelle ADU d'origine (× t) et on NE normalise PAS :
            # le master reste ainsi proportionnel à t × luminance, seule condition
            # pour que la fusion finale des masters soit juste. EXPTIME est écrit
            # dans l'en-tête → plus aucune saisie manuelle ensuite.
            exps = [float(e) for _, e in items if e and e > 0]
            uniform_t = None
            if len(exps) == len(items) and exps:
                if max(exps) <= 1.005 * min(exps):
                    uniform_t = float(np.median(exps))
            out_exptime = None
            if uniform_t is not None and mode != "mertens":
                hdr = hdr * uniform_t                    # retour à l'échelle ADU
                mx = float(np.max(hdr))
                s = mx if mx > 1.0 else 1.0              # tenir dans [0,1]
                if s > 1.0:
                    hdr = hdr / s
                out_exptime = uniform_t / s              # EXPTIME effectif
                hdr = np.maximum(hdr, 0.0).astype(np.float32)
                self.progress.emit(99,
                    f"Poses toutes à {uniform_t:g} s : échelle radiométrique "
                    f"conservée, EXPTIME={out_exptime:.6g} s écrit dans l'en-tête.")
                if abs(uniform_t - 1.0) < 1e-9:
                    self.progress.emit(99,
                        "⚠ Le temps de pose vaut 1 — s'il vient de « Estimer les "
                        "expositions » c'est une valeur RELATIVE, inutilisable pour "
                        "la fusion finale. Saisissez la vraie vitesse.")
            elif norm and mode != "mertens":
                m = float(np.percentile(hdr, 99.99))
                if m <= 0:
                    m = float(np.max(hdr))
                if m > 0:
                    v = hdr / m
                    k = float(self.p.get("knee", 0.85))
                    if k >= 0.999:
                        hdr = np.clip(v, 0.0, 1.0)
                    else:
                        hdr = np.where(
                            v <= k, v,
                            k + (1.0 - k) * np.tanh((v - k) / (1.0 - k)))
                    hdr = np.maximum(hdr, 0.0).astype(np.float32)
            self.progress.emit(100, "Fusion terminée.")
            self.done.emit(hdr.astype(np.float32), out_exptime)
        except Exception as e:
            self.error.emit(str(e))


# ── Zoom view ────────────────────────────────────────────────────────────────
class ZoomView(QGraphicsView):
    def __init__(self, scene, parent=None):
        super().__init__(scene, parent)
        self.setTransformationAnchor(QGraphicsView.ViewportAnchor.AnchorUnderMouse)
        self.setResizeAnchor(QGraphicsView.ViewportAnchor.AnchorUnderMouse)
        self.setDragMode(QGraphicsView.DragMode.ScrollHandDrag)
        self.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform)

    def wheelEvent(self, event):
        f = 1.15 if event.angleDelta().y() > 0 else 1 / 1.15
        self.scale(f, f)


# ── Main window ────────────────────────────────────────────────────────────────
class FusionHDRWindow(QMainWindow):

    def __init__(self, siril):
        super().__init__()
        self.siril = siril
        self.setWindowTitle(f"FusionHDR v1.2 (zones non saturées) — v{VERSION}")
        self.setMinimumSize(1000, 660)
        self.setStyleSheet(DARK_SS)

        # items : liste de dict {path, exptime, datamax, shape}
        self._items: list[dict] = []
        self._result: Optional[np.ndarray] = None
        self._worker: Optional[HDRWorker] = None
        self._pixitem = None
        # Dossier par défaut = dossier de travail de Siril, mémorisé/modifiable.
        self._last_dir = os.getcwd()

        self._build_ui()

    # =========================================================================
    def _build_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        h = QHBoxLayout(central)
        h.setContentsMargins(6, 6, 6, 6)
        h.setSpacing(8)

        # ── Panneau gauche ───────────────────────────────────────────────────
        left = QWidget(); left.setFixedWidth(300)
        lv = QVBoxLayout(left); lv.setContentsMargins(0, 0, 0, 0); lv.setSpacing(6)

        grp = QGroupBox("Poses (déjà alignées, linéaires)")
        gv = QVBoxLayout(grp)
        self.table = QTableWidget(0, 3)
        self.table.setHorizontalHeaderLabels(["Fichier", "Expo (s)", "Max"])
        self.table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        self.table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeMode.ResizeToContents)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setMaximumHeight(190)
        self.table.itemChanged.connect(self._on_item_edited)
        gv.addWidget(self.table)
        add_row = QHBoxLayout()
        add_row.addWidget(QLabel("Ajouter :"))
        for txt, fn in (("FITS", self._add_files), ("Séquence", self._add_sequence),
                        ("Dossier", self._add_folder)):
            b = QPushButton(txt); b.clicked.connect(fn); add_row.addWidget(b)
        add_row.addStretch()
        gv.addLayout(add_row)
        rm_row = QHBoxLayout()
        rm_row.addWidget(QLabel("Retirer :"))
        for txt, fn in (("La sélection", self._remove_sel), ("Tout", self._clear)):
            b = QPushButton(txt); b.clicked.connect(fn); rm_row.addWidget(b)
        rm_row.addStretch()
        gv.addLayout(rm_row)
        b_est = QPushButton("Estimer les expositions (recouvrement)")
        b_est.setToolTip("Calcule des rapports d'exposition relatifs à partir des\n"
                         "pixels non saturés communs — utile si les temps de pose\n"
                         "sont inconnus ou identiques.")
        b_est.clicked.connect(self._estimate_exposures)
        gv.addWidget(b_est)
        lv.addWidget(grp)

        grp2 = QGroupBox("Fusion")
        g2 = QGridLayout(grp2)
        g2.addWidget(QLabel("Mode :"), 0, 0)
        self.cmb_mode = QComboBox()
        self.cmb_mode.addItems(["Remplacement par seuil", "Mélange pondéré (SNR)",
                                "Fusion d'expositions (Mertens)",
                                "Variance minimale (MLE)",
                                "Zones non saturées (spatial)"])
        self.cmb_mode.setCurrentIndex(4)          # v1.2 : zones par défaut
        self.cmb_mode.setToolTip(
            "Zones non saturées : masque dur de saturation par pose, étendu\n"
            "d'une marge (la saturation contamine son voisinage), fondu SPATIAL\n"
            "sur la distance au masque → relais entre poses en pleine zone\n"
            "linéaire, pas de marche le long des isophotes. Pondération MLE\n"
            "parmi les poses valides.")
        g2.addWidget(self.cmb_mode, 0, 1)
        g2.addWidget(QLabel("Pleine échelle :"), 1, 0)
        self.cmb_fs = QComboBox()
        self.cmb_fs.addItems(["Auto (max des poses)", "1.0 (normalisé)", "65535 (16 bits)"])
        g2.addWidget(self.cmb_fs, 1, 1)
        g2.addWidget(QLabel("Seuil saturation :"), 2, 0)
        self.spn_sat = QDoubleSpinBox()
        self.spn_sat.setRange(0.50, 1.0); self.spn_sat.setSingleStep(0.01)
        self.spn_sat.setDecimals(2); self.spn_sat.setValue(0.85)
        g2.addWidget(self.spn_sat, 2, 1)
        g2.addWidget(QLabel("Transition (× éch.) :"), 3, 0)
        self.spn_feather = QDoubleSpinBox()
        self.spn_feather.setRange(0.0, 0.5); self.spn_feather.setSingleStep(0.02)
        self.spn_feather.setDecimals(2); self.spn_feather.setValue(0.20)
        g2.addWidget(self.spn_feather, 3, 1)
        g2.addWidget(QLabel("Marge zone saturée (px) :"), 4, 0)
        self.spn_margin = QSpinBox()
        self.spn_margin.setRange(0, 500); self.spn_margin.setValue(20)
        self.spn_margin.setToolTip("Mode Zones : exclusion étendue de N px autour de "
                                   "chaque zone saturée (bavure, halo, non-linéarité).")
        g2.addWidget(self.spn_margin, 4, 1)
        g2.addWidget(QLabel("Fondu spatial (px) :"), 5, 0)
        self.spn_sfeather = QSpinBox()
        self.spn_sfeather.setRange(1, 500); self.spn_sfeather.setValue(20)
        self.spn_sfeather.setToolTip("Mode Zones : rampe de poids sur la distance au "
                                     "masque (0 au bord de la marge → 1 à N px).")
        g2.addWidget(self.spn_sfeather, 5, 1)
        self.chk_prom = QCheckBox("Protubérances : superposition sans halo")
        self.chk_prom.setChecked(True)
        self.chk_prom.setToolTip(
            "Mode Zones : de chaque pose courte, seuls les excès localisés\n"
            "au-dessus de la ligne de base azimutale (= les protubérances,\n"
            "sans le halo) sont extraits, puis UNION (max) des N époques :\n"
            "toutes les protubérances apparaissent, aucun mélange de halos,\n"
            "plus de trous ni de demi-trous.")
        g2.addWidget(self.chk_prom, 6, 0, 1, 2)
        g2.addWidget(QLabel("   …union des N courtes :"), 7, 0)
        self.spn_promn = QSpinBox()
        self.spn_promn.setRange(1, 5); self.spn_promn.setValue(3)
        self.spn_promn.setToolTip(
            "Nombre de poses courtes dont les protubérances sont superposées.\n"
            "3 = couvre plusieurs époques (protubérances des deux bords).\n"
            "L'union (max) n'érode rien, contrairement à une médiane.")
        g2.addWidget(self.spn_promn, 7, 1)
        g2.addWidget(QLabel("   …hauteur (× rayon Lune) :"), 8, 0)
        self.spn_promh = QDoubleSpinBox()
        self.spn_promh.setRange(0.02, 1.00); self.spn_promh.setSingleStep(0.05)
        self.spn_promh.setDecimals(2); self.spn_promh.setValue(0.15)
        self.spn_promh.setToolTip(
            "Hauteur de la bande au-dessus du limbe où les protubérances sont\n"
            "cherchées, en fraction du rayon lunaire. 0,15 couvre largement les\n"
            "protubérances réelles.\n"
            "ATTENTION : au-delà de ~0,3 la bande capte les streamers de la\n"
            "couronne (eux aussi des excès azimutaux) et rehausse toute la\n"
            "couronne interne au lieu des seules protubérances.")
        g2.addWidget(self.spn_promh, 8, 1)
        self.chk_limb = QCheckBox("Limbe mono-pose (rarement utile)")
        self.chk_limb.setChecked(False)
        self.chk_limb.setToolTip(
            "Mode Zones : la bande du limbe est prise de la POSE LA PLUS COURTE\n"
            "seule. C'est une zone à règle différente : sa frontière peut se\n"
            "voir comme un cercle. Le masque lunaire continu rend cette option\n"
            "inutile dans la plupart des cas — ne l'activez qu'en dernier recours.")
        g2.addWidget(self.chk_limb, 9, 0, 1, 2)
        g2.addWidget(QLabel("Coude hautes lumières :"), 10, 0)
        self.spn_knee = QDoubleSpinBox()
        self.spn_knee.setRange(0.50, 1.0); self.spn_knee.setSingleStep(0.05)
        self.spn_knee.setDecimals(2); self.spn_knee.setValue(0.85)
        self.spn_knee.setToolTip(
            "Sortie : linéaire pur en dessous du coude, compression douce\n"
            "au-dessus (jamais d'écrêtage) → le cœur garde son modelé.\n"
            "Plus bas = cœur plus protégé. 1,00 = écrêtage dur (ancien\n"
            "comportement). Sans effet en mode Mertens.")
        g2.addWidget(self.spn_knee, 10, 1)
        lv.addWidget(grp2)

        grp3 = QGroupBox("Options")
        g3 = QGridLayout(grp3)
        self.chk_bg = QCheckBox("Soustraire le fond (ciel profond)")
        g3.addWidget(self.chk_bg, 0, 0, 1, 2)
        g3.addWidget(QLabel("Percentile fond :"), 1, 0)
        self.spn_bg = QDoubleSpinBox()
        self.spn_bg.setRange(0.0, 50.0); self.spn_bg.setValue(5.0); self.spn_bg.setDecimals(1)
        self.spn_bg.setEnabled(False)
        self.chk_bg.toggled.connect(self.spn_bg.setEnabled)
        g3.addWidget(self.spn_bg, 1, 1)
        self.chk_norm = QCheckBox("Normaliser la sortie à [0,1]")
        self.chk_norm.setChecked(True)
        g3.addWidget(self.chk_norm, 2, 0, 1, 2)
        self.chk_xcal = QCheckBox("Calibration croisée des poses (anti-bandes)")
        self.chk_xcal.setChecked(True)
        self.chk_xcal.setToolTip(
            "Ajuste gain + offset de chaque pose sur la suivante plus longue\n"
            "(régression robuste sur leurs pixels communs valides). Supprime\n"
            "les marches concentriques aux transitions entre poses (EXIF\n"
            "arrondis, non-linéarité près de la saturation, ciel variable).\n"
            "Sans effet en mode Mertens.")
        g3.addWidget(self.chk_xcal, 3, 0, 1, 2)
        lv.addWidget(grp3)

        lv.addStretch()
        self.btn_merge = QPushButton("▶  Fusionner")
        self.btn_merge.setObjectName("BtnGo")
        self.btn_merge.clicked.connect(self._merge)
        lv.addWidget(self.btn_merge)
        self.btn_save = QPushButton("Enregistrer + charger dans Siril")
        self.btn_save.setEnabled(False)
        self.btn_save.clicked.connect(self._save_load)
        lv.addWidget(self.btn_save)
        h.addWidget(left)

        # ── Panneau droit ────────────────────────────────────────────────────
        right = QWidget()
        rv = QVBoxLayout(right); rv.setContentsMargins(0, 0, 0, 0); rv.setSpacing(4)
        grp_p = QGroupBox("Aperçu HDR (autostretch)")
        pv = QVBoxLayout(grp_p)
        self._scene = QGraphicsScene()
        self._view = ZoomView(self._scene)
        pv.addWidget(self._view)
        rv.addWidget(grp_p, stretch=3)
        self.prog = QProgressBar()
        rv.addWidget(self.prog)
        self.log = QTextEdit(); self.log.setReadOnly(True); self.log.setFixedHeight(110)
        rv.addWidget(self.log)
        h.addWidget(right, stretch=1)

        self.statusBar().showMessage("Ajouter les poses (déjà alignées) à fusionner.")

    # =========================================================================
    #  Gestion des poses
    # =========================================================================
    def _add_paths(self, paths):
        added = 0
        stacked = []
        for p in paths:
            try:
                data, hdr = load_fits(p)
                exp = get_exptime(hdr) or 0.0
                n_st, lt = stack_info(hdr)
                if n_st > 1:
                    stacked.append((Path(p).name, exp, n_st, lt))
                self._items.append({
                    "path": p, "exptime": exp,
                    "datamax": float(np.max(data)), "shape": data.shape,
                })
                added += 1
            except Exception as e:
                self._log(f"✗ {Path(p).name} : {e}")
        if added:
            self._items.sort(key=lambda d: (d["exptime"], d["path"]))
            self._refresh_table()
            self._log(f"+ {added} pose(s) ajoutée(s) ({len(self._items)} au total).")
        if stacked:
            self._log(f"ℹ {len(stacked)} master(s) empilé(s) détecté(s) :")
            for nm, exp, n_st, lt in stacked:
                lts = f", LIVETIME={lt:g}" if lt else ""
                self._log(f"   {nm} : EXPTIME={exp:g}, {n_st} images empilées{lts}")
            # Siril additionne les temps de pose (cf. -nosum) alors que med/rej
            # produit une MOYENNE : le temps à utiliser ici est la pose UNITAIRE.
            susp = [s for s in stacked if s[3] and abs(s[1] - s[3]) < 1e-9 and s[2] > 1]
            if susp:
                self._log("⚠ EXPTIME = LIVETIME sur ces masters : c'est le temps de pose "
                          "CUMULÉ. Or un empilement med/rej donne une MOYENNE — saisissez "
                          "la pose UNITAIRE (EXPTIME ÷ nombre d'images), sinon les masters "
                          "seront décalés entre eux et les bandes reviendront.")

    def _add_files(self):
        files, _ = QFileDialog.getOpenFileNames(
            self, "Ajouter des poses", self._last_dir,
            "FITS (*.fit *.fits *.fts *.FIT *.FITS *.FTS)")
        if files:
            self._last_dir = str(Path(files[0]).parent)
            self._add_paths(files)

    def _add_folder(self):
        d = QFileDialog.getExistingDirectory(self, "Dossier de poses", self._last_dir)
        if not d:
            return
        self._last_dir = d
        exts = {".fit", ".fits", ".fts"}
        files = [str(f) for f in Path(d).iterdir()
                 if f.is_file() and f.suffix.lower() in exts]
        if files:
            self._add_paths(sorted(files))
        else:
            self._log("Aucun FITS dans ce dossier.")

    def _add_sequence(self):
        """Charge toutes les images d'une séquence Siril (.seq) : lit le nom de
        base puis ajoute les FITS individuels correspondants (ex. aligned_.seq)."""
        import re
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
        exts = {".fit", ".fits", ".fts"}
        files = sorted(str(p) for p in Path(f).parent.iterdir()
                       if p.is_file() and p.suffix.lower() in exts
                       and p.stem.startswith(seqname))
        if files:
            self._add_paths(files)
        else:
            self._log(f"Aucun FITS «{seqname}*» (séquence mono-fichier FITSEQ/SER ?).")

    def _remove_sel(self):
        rows = sorted({i.row() for i in self.table.selectedIndexes()}, reverse=True)
        for r in rows:
            if 0 <= r < len(self._items):
                del self._items[r]
        self._refresh_table()

    def _clear(self):
        self._items.clear()
        self._refresh_table()

    def _refresh_table(self):
        self.table.blockSignals(True)
        self.table.setRowCount(len(self._items))
        for i, it in enumerate(self._items):
            name = QTableWidgetItem(Path(it["path"]).name)
            name.setFlags(name.flags() & ~Qt.ItemFlag.ItemIsEditable)
            name.setToolTip(it["path"])
            exp = QTableWidgetItem(f"{it['exptime']:.4g}" if it["exptime"] else "")
            exp.setFlags(exp.flags() | Qt.ItemFlag.ItemIsEditable)
            mx = QTableWidgetItem(f"{it['datamax']:.4g}")
            mx.setFlags(mx.flags() & ~Qt.ItemFlag.ItemIsEditable)
            self.table.setItem(i, 0, name)
            self.table.setItem(i, 1, exp)
            self.table.setItem(i, 2, mx)
        self.table.blockSignals(False)

    def _on_item_edited(self, item):
        if item.column() != 1:
            return
        r = item.row()
        if not (0 <= r < len(self._items)):
            return
        raw = item.text()
        val = parse_exposure(raw)
        if val is None or val <= 0:
            # Ne JAMAIS revenir en arrière en silence : une saisie refusée sans
            # message laisse croire que la valeur est prise en compte, et toute
            # la fusion part alors sur de faux temps de pose.
            old = self._items[r]["exptime"]
            self.table.blockSignals(True)
            item.setText(f"{old:.6g}" if old else "")
            self.table.blockSignals(False)
            self._log(f"⚠ « {raw} » n'est pas un temps de pose valide — valeur "
                      f"inchangée ({old:.6g} s). Formats acceptés : 0.004 · 0,004 · 1/250 · 2s")
            return
        self._items[r]["exptime"] = val
        self.table.blockSignals(True)
        item.setText(f"{val:.6g}")
        self.table.blockSignals(False)
        self._log(f"Pose {Path(self._items[r]['path']).name} : {val:.6g} s")

    # =========================================================================
    #  Estimation des expositions par recouvrement
    # =========================================================================
    def _estimate_exposures(self):
        if len(self._items) < 2:
            self._log("Ajouter au moins deux poses.")
            return
        try:
            fs = self._full_scale()
            sat = self.spn_sat.value() * fs
            # Référence = pose médiane (par datamax) ; échelle ref = 1
            order = sorted(range(len(self._items)), key=lambda i: self._items[i]["datamax"])
            ref_idx = order[len(order) // 2]
            ref, _ = load_fits(self._items[ref_idx]["path"])
            ref_l = max_chan(ref)
            st = max(1, ref_l.size // 300_000)
            ref_s = ref_l.flatten()[::st]
            for i, it in enumerate(self._items):
                if i == ref_idx:
                    it["exptime"] = 1.0
                    continue
                data, _ = load_fits(it["path"])
                cur_s = max_chan(data).flatten()[::st]
                m = (ref_s > 0.02 * fs) & (ref_s < sat) & (cur_s < sat)
                if m.sum() < 50:
                    it["exptime"] = it["exptime"] or 1.0
                    continue
                k = float(np.sum(cur_s[m] * ref_s[m]) / np.sum(ref_s[m] ** 2))
                it["exptime"] = max(k, 1e-6)   # cur ≈ k·ref  →  expo relative = k
            self._items.sort(key=lambda d: (d["exptime"], d["path"]))
            self._refresh_table()
            self._log("✓ Expositions estimées (relatives, réf. = 1).")
        except Exception as e:
            self._log(f"✗ Estimation : {e}")

    # =========================================================================
    #  Fusion
    # =========================================================================
    def _full_scale(self) -> float:
        idx = self.cmb_fs.currentIndex()
        if idx == 1:
            return 1.0
        if idx == 2:
            return 65535.0
        return max((it["datamax"] for it in self._items), default=1.0) or 1.0

    def _merge(self):
        if len(self._items) < 2:
            self._log("Ajouter au moins deux poses.")
            return
        is_mertens = self.cmb_mode.currentIndex() == 2
        # Mertens est basé sur le contenu (bonne exposition) → pas besoin des
        # temps de pose ; les modes radiométriques, eux, en ont besoin.
        if not is_mertens and any(it["exptime"] <= 0 for it in self._items):
            self._log("⚠ Expositions manquantes — saisir les temps de pose ou "
                      "cliquer « Estimer les expositions » (inutile en mode Mertens).")
            return
        shapes = {it["shape"] for it in self._items}
        if len(shapes) > 1:
            self._log(f"⚠ Dimensions hétérogènes : {shapes}. Les poses doivent "
                      "être recalées au même cadrage.")
            return
        params = {
            "items":    [(it["path"], it["exptime"]) for it in self._items],
            "mode":     ("threshold", "weighted", "mertens",
                         "mle", "zones")[self.cmb_mode.currentIndex()],
            "zone_margin":  self.spn_margin.value(),
            "zone_feather": self.spn_sfeather.value(),
            "limb_mono":    self.chk_limb.isChecked(),
            "prom_overlay": self.chk_prom.isChecked(),
            "prom_n":       self.spn_promn.value(),
            "prom_h":       self.spn_promh.value(),
            "fs":       self._full_scale(),
            "sat_frac": self.spn_sat.value(),
            "feather":  self.spn_feather.value(),
            "bg_sub":   self.chk_bg.isChecked(),
            "bg_pct":   self.spn_bg.value(),
            "normalize": self.chk_norm.isChecked(),
            "knee":     self.spn_knee.value(),
            "xcal":     self.chk_xcal.isChecked(),
        }
        self.btn_merge.setEnabled(False)
        self.btn_save.setEnabled(False)
        self.prog.setValue(0)
        self._log(f"━━ Fusion ({self.cmb_mode.currentText()}, "
                  f"échelle {params['fs']:.4g}) sur {len(self._items)} poses")
        self._worker = HDRWorker(params)
        self._worker.progress.connect(lambda p, m: (self.prog.setValue(p), self._log(m)))
        self._worker.done.connect(self._on_done)
        self._worker.error.connect(lambda e: (self._log(f"✗ {e}"),
                                              self.btn_merge.setEnabled(True)))
        self._worker.start()

    def _on_done(self, result, out_exptime=None):
        self._result = result
        self._result_exptime = out_exptime
        self.btn_merge.setEnabled(True)
        self.btn_save.setEnabled(True)
        disp = autostretch_display(result)
        pm = to_pixmap(disp)
        if self._pixitem is None:
            self._pixitem = self._scene.addPixmap(pm)
        else:
            self._pixitem.setPixmap(pm)
        self._scene.setSceneRect(0, 0, pm.width(), pm.height())
        QTimer.singleShot(50, lambda: self._view.fitInView(
            self._scene.sceneRect(), Qt.AspectRatioMode.KeepAspectRatio))
        kind = "couleur" if result.ndim == 3 else "monochrome"
        self._log(f"✓ HDR {pm.width()}×{pm.height()} px ({kind}). Prêt à enregistrer.")
        self.statusBar().showMessage("Fusion terminée.")

    def _save_load(self):
        if self._result is None:
            return
        # Par défaut : le répertoire de travail de Siril (os.getcwd() = dossier
        # courant de Siril), et NON le dernier dossier des poses d'entrée — sinon
        # le HDR atterrit ailleurs et on recharge le mauvais fichier dans Corona.
        default = str(Path(os.getcwd()) / "hdr_merge.fit")
        f, _ = QFileDialog.getSaveFileName(self, "Enregistrer le HDR", default,
                                           "FITS (*.fit *.fits)")
        if not f:
            return
        try:
            save_fits(self._result, f,
                      exptime=getattr(self, '_result_exptime', None),
                      note='FusionHDR : combinaison a vitesse unique'
                           if getattr(self, '_result_exptime', None) else None)
            self._log(f"✓ Enregistré : {f}")
            stem = str(Path(f).with_suffix(""))
            try:
                # Siril en 16 bits requantifierait le HDR au chargement → la
                # couronne faible (~1e-3 de la pleine échelle) n'aurait que
                # ~65 niveaux = postérisation. On force le 32 bits d'abord.
                try:
                    self.siril.cmd("set32bits")
                    self._log("→ set32bits : Siril passé en 32 bits "
                              "(revenir avec set16bits si besoin).")
                except Exception:
                    pass
                self.siril.cmd(f'load "{stem}"')   # guillemets : chemins avec espaces
                self.siril.log(f"FusionHDR — résultat chargé : {Path(f).name}")
                self._log("✓ Chargé dans Siril.")
            except Exception:
                self._log("✓ Enregistré (chargez-le manuellement dans Siril).")
        except Exception as e:
            self._log(f"✗ {e}")

    # =========================================================================
    def _log(self, msg: str):
        self.log.append(msg)
        self.log.ensureCursorVisible()
        self.statusBar().showMessage(msg)

    def closeEvent(self, event):
        if self._worker and self._worker.isRunning():
            self._worker.abort()
            self._worker.wait(4000)
        super().closeEvent(event)


# ── Entry point ───────────────────────────────────────────────────────────────
def main():
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    siril = s.SirilInterface()
    try:
        siril.connect()
    except Exception as e:
        print(f"Avertissement : connexion Siril impossible ({e}).")
    win = FusionHDRWindow(siril)
    win.show()
    app.exec()


if __name__ == "__main__":
    main()
