##############################################
# EclipseAlign
# Alignement de brackets d'éclipse sur la Lune
# Version 1.0.0
##############################################
#
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Crédits
# -------
#   • Détection du limbe : ajustement de cercle RANSAC + Kåsa
#   • Interface : PyQt6 / conventions VeraLux

"""
EclipseAlign — recale des poses d'éclipse de luminosités très différentes
(du diamant à la couronne saturée) en les **alignant sur la Lune**.

Pourquoi
--------
L'alignement DFT de Siril corrèle le contenu global ; entre une pose très
courte (diamant, mince limbe) et une longue (couronne complète), le contenu
diffère trop → décalage faux. Le **limbe lunaire** est le seul repère commun
à toutes les poses. On l'ajuste par un cercle (RANSAC, robuste même sur un
arc partiel), puis on aligne par translation pour que les centres coïncident.

Sortie : fichiers `r_*.fit` (valeurs et en-têtes préservés, dont EXPTIME pour
FusionHDR), prêts pour la fusion HDR. Option : créer la séquence Siril.

Compatibilité
-------------
• Siril 1.3+   • Python 3.10+ (via sirilpy)
• Dépendances : numpy, astropy, scipy, PyQt6
"""

import sys
import os
import math

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
from scipy.ndimage import gaussian_filter, shift as ndshift, label, binary_closing

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
    QLabel, QPushButton, QLineEdit, QCheckBox, QSpinBox, QProgressBar,
    QTextEdit, QFileDialog, QGroupBox, QGraphicsView, QGraphicsScene,
    QGraphicsPixmapItem, QTableWidget, QTableWidgetItem, QHeaderView,
    QAbstractItemView
)
from PyQt6.QtCore import Qt, QThread, pyqtSignal, QTimer
from PyQt6.QtGui import QImage, QPixmap, QPainter, QColor, QPen

VERSION = "1.0.0"

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
QLineEdit          { background:#1e1e1e; border:1px solid #555; border-radius:3px;
                     padding:3px 6px; color:#d4d4d4; }
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
"""

# ── FITS / display helpers ───────────────────────────────────────────────────────
def load_fits(path):
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

def save_fits(data, hdr, path):
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
    g = np.ascontiguousarray(np.flipud((np.clip(lum01, 0, 1) * 255).astype(np.uint8)))
    h, w = g.shape
    return QPixmap.fromImage(QImage(g.tobytes(), w, h, w, QImage.Format.Format_Grayscale8))


# ── Détection du limbe lunaire (cercle RANSAC + refit Kåsa) ──────────────────────
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
        bright = binary_closing(bright, iterations=2)   # scelle les creux de l'anneau
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
    """Détecte la Lune. Retourne (cx, cy, r, confiance) en coords tableau
    (col, row), ou None.

    Adaptatif : si le disque sombre est **enclos** par la couronne (poses
    lumineuses/moyennes), on prend son centroïde — précis et insensible au
    halo. Sinon (diamant / poses faibles), on ajuste le **limbe par RANSAC**
    (robuste même sur un arc partiel)."""
    a = gaussian_filter(lum.astype(np.float32), 2.0)
    # Écrête le pic le plus brillant (le diamant) : sa lumière déborde dans le
    # disque et tirerait la détection du limbe vers lui.
    cap = float(np.percentile(a, 99.5))
    if cap > 0:
        a = np.minimum(a, cap)

    disk = _detect_disk_enclosed(a)
    if disk is not None:
        return disk[0], disk[1], disk[2], 1.0     # confiance maximale

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
        n = int(np.count_nonzero(d < 2.0))
        if n > best_n:
            best_n, best = n, (cx, cy, r)
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


# ── Workers ──────────────────────────────────────────────────────────────────
class DetectWorker(QThread):
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


class AlignWorker(QThread):
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
                dy = self.ref_cy - it["cy"]      # ligne
                dx = self.ref_cx - it["cx"]      # colonne
                if data.ndim == 2:
                    out = ndshift(data, (dy, dx), order=1, mode="constant", cval=0.0)
                else:
                    out = ndshift(data, (0, dy, dx), order=1, mode="constant", cval=0.0)
                # Sortie en flottant [0,1] (Siril attend du float dans [0,1]),
                # diviseur commun → échelle relative préservée pour FusionHDR.
                out = np.clip(out / self.scale, 0.0, 1.0)
                # Nommage séquence Siril : aligned_00001.fit, … (+ aligned_.seq)
                dst = str(Path(self.out_dir) / f"aligned_{n_ok+1:05d}.fit")
                save_fits(out, hdr, dst)
                n_ok += 1
                self.progress.emit(int(100*(i+1)/n),
                                   f"✓ {i+1}/{n}  —  {name}  (Δ {dx:+.1f},{dy:+.1f})")
            except Exception as e:
                self.progress.emit(int(100*(i+1)/n), f"✗ {i+1}/{n}  —  {name} : {e}")
        self.done.emit(n_ok, n, self.out_dir)


# ── Vue : clic-glissé gauche = déplacer le centre, Maj+molette = rayon, ───────────
#         molette = zoom, clic milieu glissé = déplacer la vue ────────────────────
class CenterView(QGraphicsView):
    center_at    = pyqtSignal(float, float)   # coords scène, pendant le glissé
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
        if event.modifiers() & Qt.KeyboardModifier.ShiftModifier:
            self.radius_delta.emit(1 if event.angleDelta().y() > 0 else -1)
        else:
            f = 1.15 if event.angleDelta().y() > 0 else 1 / 1.15
            self.scale(f, f)

    def mousePressEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            self._dragging = True
            p = self.mapToScene(event.position().toPoint())
            self.center_at.emit(p.x(), p.y())
            event.accept()
        elif event.button() == Qt.MouseButton.MiddleButton:
            self._pan = event.position().toPoint()
            self.setCursor(Qt.CursorShape.ClosedHandCursor); event.accept()
        else:
            super().mousePressEvent(event)

    def mouseMoveEvent(self, event):
        if self._dragging:
            p = self.mapToScene(event.position().toPoint())
            self.center_at.emit(p.x(), p.y())     # centre suivi en direct
            event.accept()
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
            self._dragging = False; event.accept()   # position validée au relâchement
        elif event.button() == Qt.MouseButton.MiddleButton:
            self._pan = None; self.setCursor(Qt.CursorShape.ArrowCursor); event.accept()
        else:
            super().mouseReleaseEvent(event)


# ── Main window ──────────────────────────────────────────────────────────────
class EclipseAlignWindow(QMainWindow):

    def __init__(self, siril):
        super().__init__()
        self.siril = siril
        self.setWindowTitle(f"EclipseAlign — v{VERSION}")
        self.setMinimumSize(1080, 680)
        self.setStyleSheet(DARK_SS)

        self._items = []        # {path, cx, cy, r, conf, shape, dmax}
        self._radius = 0.0      # rayon lunaire commun (px)
        self._cur = -1          # index affiché
        self._H = self._W = 0
        self._pixitem = None
        self._overlay = []
        # Dossier par défaut = dossier de travail de Siril, mémorisé/modifiable.
        self._last_dir = os.getcwd()
        self._detect = None
        self._align = None

        self._build_ui()

    # =========================================================================
    def _build_ui(self):
        central = QWidget(); self.setCentralWidget(central)
        h = QHBoxLayout(central); h.setContentsMargins(6, 6, 6, 6); h.setSpacing(8)

        left = QWidget(); left.setFixedWidth(330)
        lv = QVBoxLayout(left); lv.setContentsMargins(0, 0, 0, 0); lv.setSpacing(6)

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
        self.table.setMaximumHeight(220)
        self.table.itemSelectionChanged.connect(self._on_row_selected)
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
        b_det = QPushButton("↻ Détecter la Lune (toutes)")
        b_det.clicked.connect(self._detect_all)
        gv.addWidget(b_det)
        lv.addWidget(grp)

        grp2 = QGroupBox("Référence & rayon")
        g2 = QGridLayout(grp2)
        g2.addWidget(QLabel("Rayon lunaire (px) :"), 0, 0)
        self.spn_r = QSpinBox(); self.spn_r.setRange(0, 20000)
        self.spn_r.valueChanged.connect(self._on_radius_spin)
        g2.addWidget(self.spn_r, 0, 1)
        b_lock = QPushButton("Verrouiller au rayon médian")
        b_lock.clicked.connect(self._lock_median_radius)
        g2.addWidget(b_lock, 1, 0, 1, 2)
        b_all = QPushButton("⊕ Appliquer le centre affiché à toutes les poses")
        b_all.setToolTip("La Lune bouge à peine dans un bracket : réglez un centre "
                         "parfait sur une pose nette, puis propagez-le à toutes.")
        b_all.clicked.connect(self._apply_center_all)
        g2.addWidget(b_all, 2, 0, 1, 2)
        self.lbl_ref = QLabel("Référence : meilleure confiance (auto)")
        self.lbl_ref.setStyleSheet("font-size:9pt; color:#aaa;")
        self.lbl_ref.setWordWrap(True)
        g2.addWidget(self.lbl_ref, 3, 0, 1, 2)
        b_ref = QPushButton("Définir la pose affichée comme référence")
        b_ref.clicked.connect(self._set_ref_current)
        g2.addWidget(b_ref, 4, 0, 1, 2)
        lv.addWidget(grp2)

        grp3 = QGroupBox("Sortie & intégration Siril")
        g3 = QVBoxLayout(grp3)
        out_row = QHBoxLayout()
        out_row.addWidget(QLabel("Dossier :"))
        self.ed_out = QLineEdit(); self.ed_out.setPlaceholderText("vide = dossier des poses")
        out_row.addWidget(self.ed_out)
        b_o = QPushButton("…"); b_o.setFixedWidth(28); b_o.clicked.connect(self._browse_out)
        out_row.addWidget(b_o)
        g3.addLayout(out_row)
        seq_row = QHBoxLayout()
        seq_row.addWidget(QLabel("Séquence :"))
        self.chk_seq = QCheckBox("Créer aligned_.seq")
        self.chk_seq.setChecked(True)
        seq_row.addWidget(self.chk_seq)
        seq_row.addStretch()
        g3.addLayout(seq_row)
        lv.addWidget(grp3)

        hint = QLabel("Sur la pose affichée : clic-glissé = déplacer le centre "
                      "(validé au relâchement) · Maj+molette = rayon · molette = zoom · "
                      "clic milieu = déplacer la vue. Vérifiez le cercle bleu, surtout le diamant.")
        hint.setWordWrap(True); hint.setStyleSheet("font-size:8pt; color:#888;")
        lv.addWidget(hint)

        lv.addStretch()
        self._ref_idx = -1
        self.btn_align = QPushButton("▶  Aligner sur la Lune")
        self.btn_align.setObjectName("BtnGo")
        self.btn_align.clicked.connect(self._do_align)
        lv.addWidget(self.btn_align)
        h.addWidget(left)

        right = QWidget()
        rv = QVBoxLayout(right); rv.setContentsMargins(0, 0, 0, 0); rv.setSpacing(4)
        self._scene = QGraphicsScene()
        self._view = CenterView(self._scene)
        self._view.center_at.connect(self._on_center_at)
        self._view.radius_delta.connect(self._on_radius_delta)
        rv.addWidget(self._view, stretch=1)
        self.prog = QProgressBar(); rv.addWidget(self.prog)
        self.log = QTextEdit(); self.log.setReadOnly(True); self.log.setFixedHeight(96)
        rv.addWidget(self.log)
        h.addWidget(right, stretch=1)
        self.statusBar().showMessage("Ajoutez les poses du bracket.")

    # =========================================================================
    #  Liste
    # =========================================================================
    def _add_paths(self, paths):
        added = 0
        for p in paths:
            try:
                data, _ = load_fits(p)
                self._items.append({"path": p, "cx": data.shape[-1] / 2.0,
                                    "cy": data.shape[-2] / 2.0, "r": 0.0,
                                    "conf": 0.0, "shape": data.shape,
                                    "dmax": float(np.max(data))})
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
        if not d: return
        self._last_dir = d
        exts = {".fit", ".fits", ".fts"}
        files = sorted(str(f) for f in Path(d).iterdir()
                       if f.is_file() and f.suffix.lower() in exts)
        if files: self._add_paths(files)
        else: self._log("Aucun FITS dans ce dossier.")

    def _add_sequence(self):
        """Charge toutes les images d'une séquence Siril (.seq) : lit le nom de
        base puis ajoute les FITS individuels correspondants du dossier."""
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
        folder = Path(f).parent
        exts = {".fit", ".fits", ".fts"}
        files = sorted(str(p) for p in folder.iterdir()
                       if p.is_file() and p.suffix.lower() in exts
                       and p.stem.startswith(seqname))
        if files:
            self._add_paths(files)
            self._log(f"Séquence «{seqname}» : {len(files)} image(s) chargée(s).")
        else:
            self._log(f"Aucun FITS «{seqname}*» (séquence mono-fichier FITSEQ/SER ?).")

    def _remove_sel(self):
        rows = sorted({i.row() for i in self.table.selectedIndexes()}, reverse=True)
        for r in rows:
            if 0 <= r < len(self._items): del self._items[r]
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
    #  Détection
    # =========================================================================
    def _detect_all(self):
        if not self._items:
            return
        self.prog.setValue(0)
        self._log("Détection du limbe (RANSAC)…")
        paths = [it["path"] for it in self._items]
        self._detect = DetectWorker(paths)
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
        # référence auto = meilleure confiance
        if self._ref_idx < 0:
            best = max(range(len(self._items)), key=lambda k: self._items[k]["conf"],
                       default=-1)
            self._ref_idx = best
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
    #  Aperçu & édition manuelle
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
        if self._cur < 0 or self._H == 0:
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
        # Appelé en continu pendant le clic-glissé : le centre suit la souris,
        # position validée au relâchement (dernière valeur). MAJ directe des
        # cellules (pas de rebuild ni de re-sélection → pas de rechargement).
        if self._cur < 0 or self._H == 0:
            return
        it = self._items[self._cur]
        it["cx"] = float(sx)
        it["cy"] = float(self._H - 1 - sy)
        it["conf"] = 1.0                       # manuel = sûr
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
    #  Référence & alignement
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
        # Par défaut : le dossier des poses (dossier courant), pas un sous-dossier.
        out = self.ed_out.text().strip() or str(Path(self._items[0]["path"]).parent)
        self._out_dir = out
        # Diviseur commun → sortie [0,1] (corrige l'affichage « tout blanc » dans
        # Siril quand les masters sont en 16 bits) sans casser l'échelle relative.
        gmax = max((it.get("dmax", 1.0) for it in self._items), default=1.0)
        scale = gmax if gmax > 1.5 else 1.0
        self.btn_align.setEnabled(False); self.prog.setValue(0)
        self._log(f"━━ Alignement sur la Lune → {out}  "
                  f"(réf : {Path(ref['path']).name}, échelle ÷{scale:.0f})")
        self._align = AlignWorker(
            [{"path": it["path"], "cx": it["cx"], "cy": it["cy"]} for it in self._items],
            ref["cx"], ref["cy"], out, scale=scale)
        self._align.progress.connect(lambda v, m: (self.prog.setValue(v), self._log(m)))
        self._align.done.connect(self._on_align_done)
        self._align.error.connect(lambda e: (self._log(f"✗ {e}"), self.btn_align.setEnabled(True)))
        self._align.start()

    def _on_align_done(self, n_ok, n, out_dir):
        self.btn_align.setEnabled(True)
        self._log(f"✓ {n_ok}/{n} poses alignées → {out_dir}")
        try: self.siril.log(f"EclipseAlign — {n_ok}/{n} poses alignées.")
        except Exception: pass
        if self.chk_seq.isChecked() and n_ok > 0:
            self._write_seq(out_dir, n_ok)

    def _write_seq(self, out_dir, n_ok):
        """Écrit aligned_.seq directement (référence uniquement les fichiers
        aligned_NNNNN.fit, donc pas de mélange même si d'autres FITS sont
        présents dans le dossier — contrairement à `convert`)."""
        layers = 3 if (self._items and len(self._items[0]["shape"]) == 3) else 1
        ref0 = max(0, min(self._ref_idx, n_ok - 1))
        try:
            lines = ["#Siril sequence file. Généré par EclipseAlign.",
                     f"S 'aligned_' 1 {n_ok} {n_ok} 5 {ref0} 1",
                     f"L {layers}"]
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
        self.log.append(msg); self.log.ensureCursorVisible()
        self.statusBar().showMessage(msg)

    def closeEvent(self, event):
        for w in (self._detect, self._align):
            if w and w.isRunning():
                w.abort(); w.wait(4000)
        super().closeEvent(event)


def main():
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    siril = s.SirilInterface()
    try:
        siril.connect()
    except Exception as e:
        print(f"Avertissement : connexion Siril impossible ({e}).")
    win = EclipseAlignWindow(siril)
    win.show()
    app.exec()


if __name__ == "__main__":
    main()
