##############################################
# FusionHDR v1.1 — variance minimale (MLE)
# Fusion HDR radiométrique de poses bracketées
# Version 1.1.0 — variante expérimentale
##############################################
#
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Le script original FusionHDR.py reste inchangé — ceci est une variante
# à part qui ajoute le mode « Variance minimale (MLE) », mis par défaut.
#
# Crédits
# -------
#   • Interface : PyQt6 / conventions VeraLux
#   • Méthode : fusion HDR radiométrique linéaire (échelle par temps de pose)
#   • Variance minimale : estimateur du maximum de vraisemblance,
#     d'après Granados et al. 2010 (« Optimal HDR reconstruction with
#     linear digital cameras ») — poids = 1/Var(radiance), modèle de
#     bruit Var(I) = a + b·I estimé de chaque pose elle-même
#     (courbe de transfert photonique mono-image)

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
from scipy.ndimage import gaussian_filter, laplace, zoom

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QGridLayout, QLabel, QPushButton, QComboBox, QCheckBox, QSpinBox,
    QDoubleSpinBox, QProgressBar, QTextEdit, QFileDialog, QGroupBox,
    QGraphicsView, QGraphicsScene, QGraphicsPixmapItem,
    QTableWidget, QTableWidgetItem, QHeaderView, QAbstractItemView
)
from PyQt6.QtCore import Qt, QThread, pyqtSignal, QSettings, QTimer
from PyQt6.QtGui import QImage, QPixmap, QPainter

VERSION = "1.1.0"

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

def to_channels(data: np.ndarray) -> list[np.ndarray]:
    return [data] if data.ndim == 2 else [data[i] for i in range(data.shape[0])]

def max_chan(data: np.ndarray) -> np.ndarray:
    """Valeur max par pixel sur les canaux (pour le test de saturation)."""
    return data if data.ndim == 2 else data.max(axis=0)

def save_fits(data: np.ndarray, path: str):
    out = np.asarray(data, dtype=np.float32)
    fits.PrimaryHDU(out).writeto(path, overwrite=True)

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
    b, a = np.polyfit(np.asarray(bx), np.asarray(by), 1)   # pente, ordonnée
    return max(float(a), 1e-12), max(float(b), 0.0)


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
    done     = pyqtSignal(object)   # résultat float32 (H,W) ou (C,H,W)
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
            alpha, beta = np.polyfit(x[sel], y[sel], 1)
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
        A_next, B_next = 1.0, 0.0                             # référence = plus longue
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
                corr[i_short] = (A_next, B_next)
                continue
            alpha, beta = self._fit_pair(rl, rs)              # court ≈ α·long + β
            if not (0.5 < alpha < 2.0):                       # garde-fou
                alpha, beta = 1.0, 0.0
            # corrigé_court = (court − β)/α, puis correction de la pose longue
            A = A_next / alpha
            B = B_next - A_next * beta / alpha
            corr[i_short] = (A, B)
            self.progress.emit(2, f"Calibration {Path(p_s).name} ↔ {Path(p_l).name} : "
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
            if norm and mode != "mertens":
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
            self.done.emit(hdr.astype(np.float32))
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
        self.setWindowTitle(f"FusionHDR v1.1 (variance minimale) — v{VERSION}")
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
                                "Variance minimale (MLE)"])
        self.cmb_mode.setCurrentIndex(3)          # v1.1 : MLE par défaut
        self.cmb_mode.setToolTip(
            "Variance minimale : poids = 1/variance du pixel (bruit de lecture\n"
            "+ photons, mesuré sur chaque pose). Estimateur optimal, sans réglage.")
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
        g2.addWidget(QLabel("Coude hautes lumières :"), 4, 0)
        self.spn_knee = QDoubleSpinBox()
        self.spn_knee.setRange(0.50, 1.0); self.spn_knee.setSingleStep(0.05)
        self.spn_knee.setDecimals(2); self.spn_knee.setValue(0.85)
        self.spn_knee.setToolTip(
            "Sortie : linéaire pur en dessous du coude, compression douce\n"
            "au-dessus (jamais d'écrêtage) → le cœur garde son modelé.\n"
            "Plus bas = cœur plus protégé. 1,00 = écrêtage dur (ancien\n"
            "comportement). Sans effet en mode Mertens.")
        g2.addWidget(self.spn_knee, 4, 1)
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
        for p in paths:
            try:
                data, hdr = load_fits(p)
                exp = get_exptime(hdr) or 0.0
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
        try:
            self._items[r]["exptime"] = max(0.0, float(item.text().replace(",", ".")))
        except ValueError:
            self.table.blockSignals(True)
            item.setText(f"{self._items[r]['exptime']:.4g}" if self._items[r]['exptime'] else "")
            self.table.blockSignals(False)

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
                         "mle")[self.cmb_mode.currentIndex()],
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

    def _on_done(self, result):
        self._result = result
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
            save_fits(self._result, f)
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
