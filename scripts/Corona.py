##############################################
# Corona
# Rehaussement de couronne solaire (éclipse)
# Version 1.0.0
##############################################
#
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Crédits
# -------
#   • Retrait de gradient radial (normalisation azimutale)
#   • Filtre Larson–Sekanina (radial + rotationnel)
#   • Interface : PyQt6 / conventions VeraLux

"""
Corona — révèle la structure de la couronne solaire (streamers, jets,
boucles) en supprimant la forte décroissance radiale de luminosité.

Pipeline :
  1. Centre du Soleil/Lune (clic ou auto) + masque du disque lunaire
  2. Retrait du gradient radial — division par le profil azimutal moyen
     (aplatit la décroissance → la structure ressort)
  3. Filtre Larson–Sekanina — accentue les gradients radial et rotationnel
  4. Accentuation (unsharp) — détail fin
  5. Normalisation, disque noirci, chargement du résultat dans Siril

Entrée : l'image courante de Siril (typiquement le HDR produit par
FusionHDR) ou un fichier FITS. Données supposées linéaires.

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
from scipy.ndimage import gaussian_filter, gaussian_filter1d, map_coordinates

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QGridLayout, QLabel, QPushButton, QCheckBox, QSpinBox, QDoubleSpinBox,
    QProgressBar, QTextEdit, QFileDialog, QGroupBox,
    QGraphicsView, QGraphicsScene, QGraphicsPixmapItem
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
QPushButton:checked { background:#4a4a7a; border-color:#7a7aaa; }
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
QScrollBar:vertical   { background:#2b2b2b; width:10px; }
QScrollBar::handle:vertical { background:#555; border-radius:5px; }
"""

# ── Display helpers ─────────────────────────────────────────────────────────────
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

def save_fits(data, path):
    fits.PrimaryHDU(np.asarray(data, dtype=np.float32)).writeto(path, overwrite=True)


# ── Image loader (image courante de Siril) ──────────────────────────────────────
class ImageLoadWorker(QThread):
    done  = pyqtSignal(object, int, int)   # data (mono (H,W) ou (3,H,W)), W, H
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
                tmp = tempfile.mktemp(prefix="corona_")
                tmp_path = src = tmp + ".fit"
                with self.siril.image_lock():
                    self.siril.cmd(f'save "{tmp.replace(os.sep, "/")}"')
            with fits.open(src) as hdul:
                raw = hdul[0].data.astype(np.float32)
            if raw.ndim == 3 and raw.shape[2] in (3, 4) and raw.shape[0] not in (3, 4):
                raw = raw.transpose(2, 0, 1)
            if raw.ndim == 3 and raw.shape[0] == 4:
                raw = raw[:3]
            mx = float(np.max(raw))
            if mx > 1.0:
                raw = raw / (65535.0 if mx <= 65535.0 else mx)
            # mono si 3 canaux identiques
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


# ── Corona processing worker ─────────────────────────────────────────────────────
class CoronaWorker(QThread):
    progress = pyqtSignal(int, str)
    done     = pyqtSignal(object)
    error    = pyqtSignal(str)

    def __init__(self, data, params, parent=None):
        super().__init__(parent)
        self.data = data
        self.p = params

    def run(self):
        try:
            data = self.data
            p = self.p
            cx, cy = p["cx"], p["cy"]          # coords tableau (col, row)
            lunar_r = p["lunar_r"]
            H = data.shape[-2]; W = data.shape[-1]

            yy, xx = np.indices((H, W)).astype(np.float32)
            r = np.hypot(xx - cx, yy - cy)
            theta = np.arctan2(yy - cy, xx - cx)

            chans = to_channels(data)
            lum = chans[0] if len(chans) == 1 else np.mean(chans, axis=0).astype(np.float32)

            # ── 1. Modèle radial (profil azimutal moyen) ──────────────────────
            model = None
            if p["radial"]:
                self.progress.emit(15, "Retrait du gradient radial…")
                nbins = int(r.max()) + 1
                ri = np.clip(r.astype(np.int32), 0, nbins - 1)
                mask = (r >= lunar_r) & np.isfinite(lum) & (lum > 0)
                counts = np.bincount(ri[mask].ravel(), minlength=nbins).astype(np.float64)
                sums = np.bincount(ri[mask].ravel(),
                                   weights=lum[mask].ravel().astype(np.float64),
                                   minlength=nbins)
                valid = counts > 0
                prof = np.zeros(nbins)
                if valid.sum() >= 2:
                    prof[valid] = sums[valid] / counts[valid]
                    idx = np.arange(nbins)
                    prof = np.interp(idx, idx[valid], prof[valid])  # comble les trous
                    prof = gaussian_filter1d(prof, max(1.0, p["radial_smooth"]))
                model = prof[ri].astype(np.float32)
                model = np.maximum(model, 1e-6)

            def process_channel(ch):
                out = ch.astype(np.float32)
                if model is not None:
                    out = out / model
                if p["ls"]:
                    out = self._larson_sekanina(out, cx, cy, r, theta,
                                                p["ls_dr"], np.radians(p["ls_dth"]),
                                                p["ls_amount"])
                if p["unsharp"]:
                    blur = gaussian_filter(out, p["us_sigma"])
                    out = out + p["us_amount"] * (out - blur)
                return out

            results = []
            for i, ch in enumerate(chans):
                if self.isInterruptionRequested():
                    return
                self.progress.emit(40 + int(40 * i / len(chans)),
                                   f"Canal {i+1}/{len(chans)}…")
                results.append(process_channel(ch))

            self.progress.emit(90, "Normalisation…")
            stacked = results[0] if len(results) == 1 else np.stack(results)

            # Normalisation robuste sur la couronne (hors disque)
            corona = (r >= lunar_r)
            sample = (stacked[..., corona] if stacked.ndim == 3 else stacked[corona])
            lo = float(np.percentile(sample, p["clip_lo"]))
            hi = float(np.percentile(sample, p["clip_hi"]))
            if hi <= lo:
                hi = lo + 1e-6
            out = np.clip((stacked - lo) / (hi - lo), 0.0, 1.0).astype(np.float32)

            # Disque lunaire noirci
            if p["blacken"]:
                disk = r < lunar_r
                if out.ndim == 3:
                    out[:, disk] = 0.0
                else:
                    out[disk] = 0.0

            self.progress.emit(100, "Terminé.")
            self.done.emit(out)
        except Exception as e:
            self.error.emit(str(e))

    @staticmethod
    def _larson_sekanina(img, cx, cy, r, theta, dr, dth, amount):
        def samp(rr, th):
            xs = cx + rr * np.cos(th)
            ys = cy + rr * np.sin(th)
            return map_coordinates(img, [ys, xs], order=1, mode="nearest").reshape(img.shape)
        i1 = samp(r - dr, theta - dth)
        i2 = samp(r + dr, theta + dth)
        return (img + amount * (2.0 * img - i1 - i2)).astype(np.float32)


# ── Scene with center pick ───────────────────────────────────────────────────────
class PickScene(QGraphicsScene):
    clicked = pyqtSignal(float, float)

    def mousePressEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            p = event.scenePos()
            self.clicked.emit(p.x(), p.y())
        super().mousePressEvent(event)


class PickView(QGraphicsView):
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
class CoronaWindow(QMainWindow):

    def __init__(self, siril):
        super().__init__()
        self.siril = siril
        self.setWindowTitle(f"Corona — v{VERSION}")
        self.setMinimumSize(1020, 680)
        self.setStyleSheet(DARK_SS)

        self.data = None            # (H,W) ou (3,H,W) linéaire
        self.W = self.H = 0
        self.cx = self.cy = 0.0      # coords tableau (col, row depuis le bas)
        self.result = None
        self._loader = None
        self._worker = None
        self._pixitem = None
        self._overlay = []

        self._build_ui()
        self._load_image()

    # =========================================================================
    def _build_ui(self):
        central = QWidget(); self.setCentralWidget(central)
        h = QHBoxLayout(central); h.setContentsMargins(6, 6, 6, 6); h.setSpacing(8)

        left = QWidget(); left.setFixedWidth(260)
        lv = QVBoxLayout(left); lv.setContentsMargins(0, 0, 0, 0); lv.setSpacing(6)

        # Source
        row0 = QHBoxLayout()
        b_reload = QPushButton("⟳ Image Siril"); b_reload.clicked.connect(self._load_image)
        b_open = QPushButton("Fichier…"); b_open.clicked.connect(self._open_file)
        row0.addWidget(b_reload); row0.addWidget(b_open)
        lv.addLayout(row0)

        # Centre
        grp_c = QGroupBox("Centre & disque")
        gc = QGridLayout(grp_c)
        self.btn_pick = QPushButton("◎ Cliquer le centre")
        self.btn_pick.setCheckable(True)
        self.btn_pick.toggled.connect(self._toggle_pick)
        gc.addWidget(self.btn_pick, 0, 0, 1, 2)
        b_auto = QPushButton("Auto"); b_auto.clicked.connect(self._auto_center)
        gc.addWidget(b_auto, 1, 0)
        self.lbl_center = QLabel("centre : —")
        self.lbl_center.setStyleSheet("font-size:9pt; color:#aaa;")
        gc.addWidget(self.lbl_center, 1, 1)
        gc.addWidget(QLabel("Rayon lunaire (px) :"), 2, 0)
        self.spn_lunar = QSpinBox(); self.spn_lunar.setRange(0, 20000); self.spn_lunar.setValue(80)
        self.spn_lunar.valueChanged.connect(self._draw_overlay)
        gc.addWidget(self.spn_lunar, 2, 1)
        lv.addWidget(grp_c)

        # Gradient radial
        grp_r = QGroupBox("Retrait du gradient radial")
        gr = QGridLayout(grp_r)
        self.chk_radial = QCheckBox("Activé"); self.chk_radial.setChecked(True)
        gr.addWidget(self.chk_radial, 0, 0, 1, 2)
        gr.addWidget(QLabel("Lissage profil :"), 1, 0)
        self.spn_rsmooth = QDoubleSpinBox(); self.spn_rsmooth.setRange(1, 200)
        self.spn_rsmooth.setValue(8.0)
        gr.addWidget(self.spn_rsmooth, 1, 1)
        lv.addWidget(grp_r)

        # Larson-Sekanina
        grp_ls = QGroupBox("Larson–Sekanina")
        gl = QGridLayout(grp_ls)
        self.chk_ls = QCheckBox("Activé")
        gl.addWidget(self.chk_ls, 0, 0, 1, 2)
        gl.addWidget(QLabel("Δ radial (px) :"), 1, 0)
        self.spn_dr = QDoubleSpinBox(); self.spn_dr.setRange(0, 100); self.spn_dr.setValue(0.0)
        gl.addWidget(self.spn_dr, 1, 1)
        gl.addWidget(QLabel("Δ angle (°) :"), 2, 0)
        self.spn_dth = QDoubleSpinBox(); self.spn_dth.setRange(0, 45); self.spn_dth.setValue(5.0)
        gl.addWidget(self.spn_dth, 2, 1)
        gl.addWidget(QLabel("Intensité :"), 3, 0)
        self.spn_lsamt = QDoubleSpinBox(); self.spn_lsamt.setRange(0, 5); self.spn_lsamt.setValue(1.0)
        self.spn_lsamt.setSingleStep(0.1)
        gl.addWidget(self.spn_lsamt, 3, 1)
        lv.addWidget(grp_ls)

        # Unsharp
        grp_u = QGroupBox("Accentuation (unsharp)")
        gu = QGridLayout(grp_u)
        self.chk_us = QCheckBox("Activé")
        gu.addWidget(self.chk_us, 0, 0, 1, 2)
        gu.addWidget(QLabel("Sigma (px) :"), 1, 0)
        self.spn_ussig = QDoubleSpinBox(); self.spn_ussig.setRange(0.5, 100); self.spn_ussig.setValue(3.0)
        gu.addWidget(self.spn_ussig, 1, 1)
        gu.addWidget(QLabel("Intensité :"), 2, 0)
        self.spn_usamt = QDoubleSpinBox(); self.spn_usamt.setRange(0, 5); self.spn_usamt.setValue(1.0)
        self.spn_usamt.setSingleStep(0.1)
        gu.addWidget(self.spn_usamt, 2, 1)
        lv.addWidget(grp_u)

        self.chk_black = QCheckBox("Noircir le disque lunaire"); self.chk_black.setChecked(True)
        lv.addWidget(self.chk_black)

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
        self._scene = PickScene()
        self._scene.clicked.connect(self._on_scene_click)
        self._view = PickView(self._scene)
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
        f, _ = QFileDialog.getOpenFileName(self, "Ouvrir un FITS", str(Path.home()),
                                           "FITS (*.fit *.fits *.fts *.FIT *.FITS *.FTS)")
        if f:
            self._load_image(f)

    def _on_image(self, data, W, H):
        self.data = data
        self.W, self.H = W, H
        self.result = None
        self.btn_save.setEnabled(False)
        self._show(autostretch_display(data))
        if self.cx == 0 and self.cy == 0:
            self._auto_center()
        else:
            self._draw_overlay()
        self._log(f"Image {W}×{H} ({'couleur' if data.ndim == 3 else 'mono'}) chargée.")

    def _show(self, disp):
        pm = to_pixmap(disp)
        if self._pixitem is None:
            self._pixitem = self._scene.addPixmap(pm)
        else:
            self._pixitem.setPixmap(pm)
        self._scene.setSceneRect(0, 0, pm.width(), pm.height())
        QTimer.singleShot(50, lambda: self._view.fitInView(
            self._scene.sceneRect(), Qt.AspectRatioMode.KeepAspectRatio))

    # =========================================================================
    #  Centre
    # =========================================================================
    def _toggle_pick(self, on):
        self._view.setDragMode(QGraphicsView.DragMode.NoDrag if on
                               else QGraphicsView.DragMode.ScrollHandDrag)
        self.statusBar().showMessage("Cliquez le centre du Soleil." if on else "")

    def _on_scene_click(self, sx, sy):
        if not self.btn_pick.isChecked() or self.data is None:
            return
        self.cx = float(sx)
        self.cy = float(self.H - 1 - sy)     # scène → tableau (depuis le bas)
        self.btn_pick.setChecked(False)
        self._update_center_label()
        self._draw_overlay()

    def _auto_center(self):
        if self.data is None:
            return
        lum = self.data if self.data.ndim == 2 else np.mean(self.data, axis=0)
        thr = float(np.median(lum))
        w = np.clip(lum - thr, 0, None)
        tot = float(w.sum())
        if tot <= 0:
            self.cx, self.cy = self.W / 2.0, self.H / 2.0
        else:
            yy, xx = np.indices(lum.shape)
            self.cx = float((w * xx).sum() / tot)
            self.cy = float((w * yy).sum() / tot)
        self._update_center_label()
        self._draw_overlay()
        self._log(f"Centre auto : ({self.cx:.0f}, {self.cy:.0f})")

    def _update_center_label(self):
        self.lbl_center.setText(f"centre : ({self.cx:.0f}, {self.cy:.0f})")

    def _draw_overlay(self):
        for it in self._overlay:
            self._scene.removeItem(it)
        self._overlay.clear()
        if self.data is None:
            return
        sx = self.cx
        sy = self.H - 1 - self.cy            # tableau → scène
        pen = QPen(QColor("#ff3333")); pen.setCosmetic(True); pen.setWidth(2)
        L = 14
        self._overlay.append(self._scene.addLine(sx - L, sy, sx + L, sy, pen))
        self._overlay.append(self._scene.addLine(sx, sy - L, sx, sy + L, pen))
        lr = self.spn_lunar.value()
        if lr > 0:
            penc = QPen(QColor("#33aaff")); penc.setCosmetic(True); penc.setWidth(2)
            self._overlay.append(self._scene.addEllipse(sx - lr, sy - lr, 2 * lr, 2 * lr, penc))

    # =========================================================================
    #  Traitement
    # =========================================================================
    def _apply(self):
        if self.data is None:
            self._log("Aucune image.")
            return
        params = {
            "cx": self.cx, "cy": self.cy,
            "lunar_r": float(self.spn_lunar.value()),
            "radial": self.chk_radial.isChecked(),
            "radial_smooth": self.spn_rsmooth.value(),
            "ls": self.chk_ls.isChecked(),
            "ls_dr": self.spn_dr.value(),
            "ls_dth": self.spn_dth.value(),
            "ls_amount": self.spn_lsamt.value(),
            "unsharp": self.chk_us.isChecked(),
            "us_sigma": self.spn_ussig.value(),
            "us_amount": self.spn_usamt.value(),
            "clip_lo": 0.5, "clip_hi": 99.7,
            "blacken": self.chk_black.isChecked(),
        }
        self.btn_apply.setEnabled(False)
        self.btn_save.setEnabled(False)
        self.prog.setValue(0)
        self._log("━━ Traitement de la couronne…")
        self._worker = CoronaWorker(self.data, params)
        self._worker.progress.connect(lambda p, m: (self.prog.setValue(p), self._log(m)))
        self._worker.done.connect(self._on_done)
        self._worker.error.connect(lambda e: (self._log(f"✗ {e}"),
                                              self.btn_apply.setEnabled(True)))
        self._worker.start()

    def _on_done(self, result):
        self.result = result
        self.btn_apply.setEnabled(True)
        self.btn_save.setEnabled(True)
        self._show(autostretch_display(result))
        self._draw_overlay()
        self._log("✓ Couronne traitée. Ajustez les paramètres ou enregistrez.")

    def _save_load(self):
        if self.result is None:
            return
        default = str(Path.home() / "corona.fit")
        f, _ = QFileDialog.getSaveFileName(self, "Enregistrer", default, "FITS (*.fit *.fits)")
        if not f:
            return
        try:
            save_fits(self.result, f)
            self._log(f"✓ Enregistré : {f}")
            try:
                self.siril.cmd("load", str(Path(f).with_suffix("")))
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
