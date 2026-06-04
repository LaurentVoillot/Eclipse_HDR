##############################################
# FusionHDR
# Fusion HDR radiométrique de poses bracketées
# Version 1.0.0
##############################################
#
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Crédits
# -------
#   • Interface : PyQt6 / conventions VeraLux
#   • Méthode : fusion HDR radiométrique linéaire (échelle par temps de pose)

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
  • Remplacement par seuil (défaut) — la plus longue pose non saturée gagne
  • Mélange pondéré — moyenne pondérée par l'exposition (SNR), hors saturation

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

s.ensure_installed("numpy", "astropy", "PyQt6")

import numpy as np
from pathlib import Path
from typing import Optional

from astropy.io import fits

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QGridLayout, QLabel, QPushButton, QComboBox, QCheckBox, QSpinBox,
    QDoubleSpinBox, QProgressBar, QTextEdit, QFileDialog, QGroupBox,
    QGraphicsView, QGraphicsScene, QGraphicsPixmapItem,
    QTableWidget, QTableWidgetItem, QHeaderView, QAbstractItemView
)
from PyQt6.QtCore import Qt, QThread, pyqtSignal, QSettings, QTimer
from PyQt6.QtGui import QImage, QPixmap, QPainter

VERSION = "1.0.0"

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

    def run(self):
        try:
            items   = self.p["items"]          # [(path, exptime), …]
            mode    = self.p["mode"]            # "threshold" | "weighted"
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

            if mode == "threshold":
                # Tri par exposition croissante : on part de la plus courte
                # (la moins saturée) et chaque pose plus longue remplace là où
                # elle n'est pas saturée → la plus longue non saturée gagne.
                order = sorted(range(n), key=lambda i: items[i][1])
                p0, e0 = items[order[0]]
                self.progress.emit(int(100*1/(n+1)), f"Base : {Path(p0).name}")
                base, _ = load_fits(p0)
                hdr = self._radiance(base, e0, bg_sub, bg_pct)
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
                    if hdr.ndim == 3:
                        w = w[None, :, :]
                    hdr = w * rad + (1.0 - w) * hdr
                    self.progress.emit(int(100*(k)/(n+1)), f"{k}/{n}  —  {Path(path).name}")

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

            # Normalisation de sortie
            if norm:
                m = float(np.max(hdr))
                if m > 0:
                    hdr = hdr / m
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
        self.setWindowTitle(f"FusionHDR — v{VERSION}")
        self.setMinimumSize(1000, 660)
        self.setStyleSheet(DARK_SS)

        # items : liste de dict {path, exptime, datamax, shape}
        self._items: list[dict] = []
        self._result: Optional[np.ndarray] = None
        self._worker: Optional[HDRWorker] = None
        self._pixitem = None

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
        row = QHBoxLayout()
        b_add = QPushButton("Ajouter…");  b_add.clicked.connect(self._add_files)
        b_fld = QPushButton("Dossier…");  b_fld.clicked.connect(self._add_folder)
        b_del = QPushButton("Retirer");   b_del.clicked.connect(self._remove_sel)
        b_clr = QPushButton("Vider");     b_clr.clicked.connect(self._clear)
        for b in (b_add, b_fld, b_del, b_clr):
            row.addWidget(b)
        gv.addLayout(row)
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
        self.cmb_mode.addItems(["Remplacement par seuil", "Mélange pondéré (SNR)"])
        g2.addWidget(self.cmb_mode, 0, 1)
        g2.addWidget(QLabel("Pleine échelle :"), 1, 0)
        self.cmb_fs = QComboBox()
        self.cmb_fs.addItems(["Auto (max des poses)", "1.0 (normalisé)", "65535 (16 bits)"])
        g2.addWidget(self.cmb_fs, 1, 1)
        g2.addWidget(QLabel("Seuil saturation :"), 2, 0)
        self.spn_sat = QDoubleSpinBox()
        self.spn_sat.setRange(0.50, 1.0); self.spn_sat.setSingleStep(0.01)
        self.spn_sat.setDecimals(2); self.spn_sat.setValue(0.95)
        g2.addWidget(self.spn_sat, 2, 1)
        g2.addWidget(QLabel("Transition (× éch.) :"), 3, 0)
        self.spn_feather = QDoubleSpinBox()
        self.spn_feather.setRange(0.0, 0.5); self.spn_feather.setSingleStep(0.02)
        self.spn_feather.setDecimals(2); self.spn_feather.setValue(0.10)
        g2.addWidget(self.spn_feather, 3, 1)
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
            self, "Ajouter des poses", str(Path.home()),
            "FITS (*.fit *.fits *.fts *.FIT *.FITS *.FTS)")
        if files:
            self._add_paths(files)

    def _add_folder(self):
        d = QFileDialog.getExistingDirectory(self, "Dossier de poses", str(Path.home()))
        if not d:
            return
        exts = {".fit", ".fits", ".fts"}
        files = [str(f) for f in Path(d).iterdir()
                 if f.is_file() and f.suffix.lower() in exts]
        if files:
            self._add_paths(sorted(files))
        else:
            self._log("Aucun FITS dans ce dossier.")

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
        if any(it["exptime"] <= 0 for it in self._items):
            self._log("⚠ Expositions manquantes — saisir les temps de pose ou "
                      "cliquer « Estimer les expositions ».")
            return
        shapes = {it["shape"] for it in self._items}
        if len(shapes) > 1:
            self._log(f"⚠ Dimensions hétérogènes : {shapes}. Les poses doivent "
                      "être recalées au même cadrage.")
            return
        params = {
            "items":    [(it["path"], it["exptime"]) for it in self._items],
            "mode":     "threshold" if self.cmb_mode.currentIndex() == 0 else "weighted",
            "fs":       self._full_scale(),
            "sat_frac": self.spn_sat.value(),
            "feather":  self.spn_feather.value(),
            "bg_sub":   self.chk_bg.isChecked(),
            "bg_pct":   self.spn_bg.value(),
            "normalize": self.chk_norm.isChecked(),
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
        default = str(Path(self._items[0]["path"]).parent / "hdr_merge.fit")
        f, _ = QFileDialog.getSaveFileName(self, "Enregistrer le HDR", default,
                                           "FITS (*.fit *.fits)")
        if not f:
            return
        try:
            save_fits(self._result, f)
            self._log(f"✓ Enregistré : {f}")
            stem = str(Path(f).with_suffix(""))
            try:
                self.siril.cmd("load", stem)
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
