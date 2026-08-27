##############################################
# ConvertParVitesse
# Conversion RAW → FITS, rangés par vitesse d'obturation
# Version 1.0.0
##############################################
#
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Crédits
# -------
#   • Lecture EXIF : exifread (repli sur exiftool si présent)
#   • Conversion / dématriçage : Siril (commande `convert -debayer`)
#   • Interface : PyQt6 / conventions VeraLux

"""
ConvertParVitesse — prépare une session bracketée pour l'alignement.

Cas d'usage : totalité d'éclipse (ou diamants, ou toute série où seule la
VITESSE varie — ouverture et ISO constants). Chaque vitesse doit être
alignée et empilée séparément avant la fusion HDR : il faut donc d'abord
séparer les poses par temps d'exposition.

Ce que fait le script
---------------------
  1. lit le temps de pose EXIF de chaque RAW (sans rien convertir) ;
  2. regroupe les poses par vitesse (tolérance réglable) et affiche le
     plan de rangement — cochez/décochez les groupes à traiter (utile
     pour ne prendre qu'une partie, ex. les diamants) ;
  3. crée un sous-dossier par vitesse, nommé comme la vitesse avec « / »
     remplacé par « - » (1/160 s → « 1-160 », 2 s → « 2s ») ;
  4. y convertit les RAW en FITS **avec dématriçage**, via Siril, ce qui
     crée aussi la séquence `.seq` prête pour l'alignement.

Arborescence obtenue :

    totalité/
      raw/            ← vos RAW (inchangés)
      1-4000/         ← FITS + séquence de la vitesse 1/4000 s
      1-160/
      2s/

Les RAW ne sont ni déplacés ni copiés : des liens temporaires sont posés
dans le dossier cible le temps de la conversion, puis retirés.

Compatibilité
-------------
• Siril 1.3+   • Python 3.10+ (via sirilpy)
• Dépendances : exifread, PyQt6
• Formats : tout ce que Siril/libraw sait lire. La lecture EXIF couvre
  les RAW TIFF (CR2, NEF, ARW, DNG, ORF, PEF, RAF…) ; pour les
  conteneurs récents (CR3…), installez `exiftool` — le script s'en sert
  automatiquement s'il est présent.
"""

import sys
import os
import html
import shutil
import subprocess
from fractions import Fraction

try:
    import sirilpy as s
except ImportError:
    print("Erreur : module sirilpy introuvable.")
    sys.exit(1)

s.ensure_installed("exifread", "PyQt6")

from pathlib import Path
from typing import Optional

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
    QLabel, QLineEdit, QPushButton, QComboBox, QCheckBox,
    QDoubleSpinBox, QProgressBar, QTextEdit, QFileDialog, QGroupBox,
    QTableWidget, QTableWidgetItem, QHeaderView, QAbstractItemView
)
from PyQt6.QtCore import Qt, QThread, pyqtSignal, QSettings
from PyQt6.QtGui import QColor

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
QPushButton#BtnScan { background:#1e3d6b; color:#fff; border:1px solid #2a60a0;
                      padding:5px 10px; }
QPushButton#BtnScan:hover    { background:#255080; }
QPushButton#BtnScan:disabled { background:#333; color:#666; border-color:#444; }
QLineEdit          { background:#1e1e1e; border:1px solid #555; border-radius:3px;
                     padding:3px 6px; color:#d4d4d4; }
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
QTableWidget       { background:#1e1e1e; gridline-color:#3a3a3a; color:#d4d4d4; }
QTableWidget::item:selected { background:#3a5f3a; }
QHeaderView::section { background:#333; color:#aaa; padding:3px; border:1px solid #444; }
QScrollBar:vertical   { background:#2b2b2b; width:10px; }
QScrollBar::handle:vertical { background:#555; border-radius:5px; }
"""

# Extensions RAW / images que Siril (libraw) sait convertir.
RAW_EXT = {
    ".cr2", ".cr3", ".crw", ".nef", ".nrw", ".arw", ".srf", ".sr2", ".dng",
    ".orf", ".pef", ".ptx", ".raf", ".rw2", ".raw", ".rwl", ".mrw", ".srw",
    ".x3f", ".3fr", ".fff", ".iiq", ".mos", ".kdc", ".dcr", ".erf",
    ".tif", ".tiff", ".jpg", ".jpeg", ".png", ".bmp",
}


# ── Lecture du temps de pose ──────────────────────────────────────────────────
def _exposure_exifread(path: str):
    """→ (num, den, secondes) ou None."""
    try:
        import exifread
    except Exception:
        return None
    try:
        with open(path, "rb") as fh:
            tags = exifread.process_file(fh, details=False)
    except Exception:
        return None
    for key in ("EXIF ExposureTime", "Image ExposureTime"):
        t = tags.get(key)
        if t is not None and getattr(t, "values", None):
            try:
                r = t.values[0]
                num, den = int(r.num), int(r.den)
                if num > 0 and den > 0:
                    return num, den, num / den
            except Exception:
                pass
    return None


def _exposure_exiftool(path: str, exe: str):
    """Repli pour les conteneurs non-TIFF (CR3…) → (num, den, secondes)."""
    try:
        out = subprocess.run([exe, "-s3", "-n", "-ExposureTime", path],
                             capture_output=True, text=True, timeout=30)
        v = out.stdout.strip()
        if not v:
            return None
        sec = float(v)
        if sec <= 0:
            return None
        fr = Fraction(sec).limit_denominator(100000)
        return fr.numerator, fr.denominator, sec
    except Exception:
        return None


def speed_label(sec: float, num: Optional[int] = None,
                den: Optional[int] = None) -> str:
    """Nom de dossier pour une vitesse : 1/160 s → « 1-160 », 2 s → « 2s »,
    1,3 s → « 1s3 ». Jamais de « / » ni de « . » (noms de fichiers portables)."""
    if sec >= 1.0:
        if abs(sec - round(sec)) < 1e-6:
            return f"{int(round(sec))}s"
        txt = f"{sec:.2f}".rstrip("0").rstrip(".")
        return txt.replace(".", "s")
    if num == 1 and den and den > 1:
        return f"1-{den}"
    inv = 1.0 / sec
    return f"1-{int(round(inv))}" if inv >= 1 else f"1-{inv:.3g}".replace(".", "s")


def group_by_speed(entries, tol_pct: float):
    """Regroupe des poses triées par temps croissant.

    Deux poses sont dans le même groupe si leur écart relatif est sous la
    tolérance (par défaut 3 % : fusionne les arrondis APEX du type 1/1000
    vs 1/1024, sans mélanger deux crans de 1/3 de diaph, distants de 26 %).
    Retourne [{label, sec, num, den, files:[...]}] + la liste des inconnus.
    """
    known = [e for e in entries if e.get("sec")]
    unknown = [e for e in entries if not e.get("sec")]
    known.sort(key=lambda e: (e["sec"], e["path"]))
    groups = []
    for e in known:
        if groups:
            ref = groups[-1]["sec_ref"]
            if abs(e["sec"] - ref) <= (tol_pct / 100.0) * ref:
                groups[-1]["files"].append(e)
                continue
        groups.append({"sec_ref": e["sec"], "files": [e]})
    out = []
    for g in groups:
        mid = g["files"][len(g["files"]) // 2]          # pose médiane du groupe
        out.append({
            "label": speed_label(mid["sec"], mid.get("num"), mid.get("den")),
            "sec": mid["sec"], "num": mid.get("num"), "den": mid.get("den"),
            "files": g["files"],
        })
    return out, unknown


def speed_text(sec: float, num: Optional[int], den: Optional[int]) -> str:
    """Libellé lisible : « 1/160 s », « 2 s »."""
    if sec >= 1.0:
        return f"{sec:g} s"
    if num == 1 and den:
        return f"1/{den} s"
    return f"1/{round(1.0/sec)} s"


def place_file(src: str, dst: Path, mode: str) -> str:
    """Met le RAW à disposition dans le dossier cible : lien symbolique, à
    défaut lien dur, à défaut copie. Retourne la méthode réellement utilisée."""
    if mode == "copie":
        shutil.copy2(src, dst)
        return "copie"
    try:
        os.symlink(os.path.abspath(src), dst)
        return "lien"
    except Exception:
        pass
    try:
        os.link(src, dst)
        return "lien dur"
    except Exception:
        pass
    shutil.copy2(src, dst)
    return "copie"


def siril_safe_cmd(siril, *args) -> bool:
    try:
        siril.cmd(*args)
        return True
    except Exception:
        return False


def siril_safe_log(siril, msg: str):
    try:
        siril.log(msg)
    except Exception:
        pass


# ── Workers ───────────────────────────────────────────────────────────────────
class ScanWorker(QThread):
    progress = pyqtSignal(int, str)
    done     = pyqtSignal(object)
    error    = pyqtSignal(str)

    def __init__(self, files, parent=None):
        super().__init__(parent)
        self.files = files
        self._abort = False

    def abort(self):
        self._abort = True

    def run(self):
        exe = shutil.which("exiftool")
        n = len(self.files)
        entries = []
        n_tool = 0
        for i, p in enumerate(self.files):
            if self._abort:
                return
            r = _exposure_exifread(p)
            if r is None and exe:
                r = _exposure_exiftool(p, exe)
                if r is not None:
                    n_tool += 1
            entries.append({"path": p,
                            "num": r[0] if r else None,
                            "den": r[1] if r else None,
                            "sec": r[2] if r else None})
            if (i + 1) % 10 == 0 or i + 1 == n:
                self.progress.emit(int(100 * (i + 1) / n),
                                   f"Lecture EXIF {i+1}/{n}  —  {Path(p).name}")
        if n_tool:
            self.progress.emit(100, f"({n_tool} fichier(s) lus via exiftool)")
        self.done.emit(entries)


class ConvertWorker(QThread):
    progress = pyqtSignal(int, str)
    done     = pyqtSignal(int, int, bool)      # groupes OK, total, interrompu
    error    = pyqtSignal(str)

    def __init__(self, siril, params, parent=None):
        super().__init__(parent)
        self.siril = siril
        self.p = params
        self._abort = False

    def abort(self):
        self._abort = True

    def run(self):
        p = self.p
        groups = p["groups"]
        out_root = Path(p["out_root"])
        prefix = p["prefix"]
        mode = p["place_mode"]
        debayer = p["debayer"]
        purge_all = p["purge_all"]
        orig_wd = os.getcwd()
        n = len(groups)
        n_ok = 0

        if p["set32"]:
            if siril_safe_cmd(self.siril, "set32bits"):
                self.progress.emit(0, "→ Siril passé en 32 bits.")

        try:
            for gi, g in enumerate(groups):
                if self._abort:
                    break
                label = g["label"]
                base = f"{prefix}{label}" if prefix else label
                target = out_root / label
                pc0 = int(100 * gi / n)
                self.progress.emit(pc0, f"━━ {speed_text(g['sec'], g.get('num'), g.get('den'))} "
                                        f"→ {target.name}/  ({len(g['files'])} images)")
                try:
                    target.mkdir(parents=True, exist_ok=True)
                except Exception as e:
                    self.progress.emit(pc0, f"✗ {target} : {e}")
                    continue

                # Purge : nos propres sorties d'un run précédent (sinon `convert`
                # ré-intégrerait ces FITS dans la nouvelle séquence).
                try:
                    for pat in (f"{base}_*.fit", f"{base}_*.fits", f"{base}*.seq"):
                        for old in target.glob(pat):
                            old.unlink()
                    if purge_all:
                        for pat in ("*.fit", "*.fits", "*.seq"):
                            for old in target.glob(pat):
                                old.unlink()
                except Exception as e:
                    self.progress.emit(pc0, f"⚠ Purge partielle : {e}")

                # FITS étrangers restants → `convert` les avalerait : on saute.
                strays = [f for f in target.iterdir()
                          if f.is_file() and f.suffix.lower() in (".fit", ".fits")]
                if strays:
                    self.progress.emit(pc0,
                        f"⚠ {len(strays)} FITS déjà présents dans {target.name}/ et non "
                        "produits par ce script — groupe ignoré. Cochez « Vider le "
                        "dossier cible » ou choisissez un autre dossier de sortie.")
                    continue

                # Liens temporaires, numérotés dans l'ordre chronologique.
                links = []
                used = set()
                try:
                    for i, f in enumerate(g["files"]):
                        if self._abort:
                            break
                        ext = Path(f["path"]).suffix.lower()
                        dst = target / f"_src_{i+1:05d}{ext}"
                        used.add(place_file(f["path"], dst, mode))
                        links.append(dst)
                except Exception as e:
                    self.progress.emit(pc0, f"✗ Mise en place : {e}")
                    for l in links:
                        try: l.unlink()
                        except Exception: pass
                    continue

                if self._abort:
                    for l in links:
                        try: l.unlink()
                        except Exception: pass
                    break

                # Conversion par Siril, dans le dossier cible.
                ok = siril_safe_cmd(self.siril, f'cd "{target.resolve()}"')
                if not ok:
                    self.progress.emit(pc0, f"✗ cd impossible : {target}")
                else:
                    args = ["convert", base] + (["-debayer"] if debayer else [])
                    ok = siril_safe_cmd(self.siril, *args)
                    if not ok and debayer:
                        self.progress.emit(pc0, "⚠ convert -debayer a échoué — "
                                                "nouvel essai sans dématriçage.")
                        ok = siril_safe_cmd(self.siril, "convert", base)

                # Retrait des liens : ne restent que les FITS + la séquence.
                for l in links:
                    try: l.unlink()
                    except Exception: pass

                made = sorted(target.glob(f"{base}_*.fit")) + \
                       sorted(target.glob(f"{base}_*.fits"))
                if ok and made:
                    n_ok += 1
                    via = "/".join(sorted(used)) or "—"
                    self.progress.emit(int(100 * (gi + 1) / n),
                        f"✓ {target.name}/ : {len(made)} FITS + séquence «{base}» ({via}).")
                else:
                    self.progress.emit(int(100 * (gi + 1) / n),
                        f"✗ {target.name}/ : conversion sans résultat.")
        finally:
            siril_safe_cmd(self.siril, f'cd "{orig_wd}"')

        self.done.emit(n_ok, n, self._abort)


# ── Fenêtre principale ────────────────────────────────────────────────────────
class ConvertWindow(QMainWindow):

    def __init__(self, siril):
        super().__init__()
        self.siril = siril
        self.setWindowTitle(f"Convert par vitesse — v{VERSION}")
        self.setMinimumSize(980, 620)
        self.setStyleSheet(DARK_SS)

        self._entries = []
        self._groups = []
        self._unknown = []
        self._scan = None
        self._conv = None

        self._build_ui()
        self._restore_settings()
        self._guess_folders()

    # =========================================================================
    def _build_ui(self):
        central = QWidget(); self.setCentralWidget(central)
        h = QHBoxLayout(central); h.setContentsMargins(6, 6, 6, 6); h.setSpacing(8)

        left = QWidget(); left.setFixedWidth(340)
        lv = QVBoxLayout(left); lv.setContentsMargins(0, 0, 0, 0); lv.setSpacing(6)

        grp_d = QGroupBox("Dossiers")
        gd = QGridLayout(grp_d)
        gd.addWidget(QLabel("RAW :"), 0, 0)
        self.ed_raw = QLineEdit()
        self.ed_raw.setPlaceholderText("dossier contenant les RAW")
        gd.addWidget(self.ed_raw, 0, 1)
        b1 = QPushButton("…"); b1.setFixedWidth(28); b1.clicked.connect(self._browse_raw)
        gd.addWidget(b1, 0, 2)
        gd.addWidget(QLabel("Sortie :"), 1, 0)
        self.ed_out = QLineEdit()
        self.ed_out.setPlaceholderText("dossier où créer les sous-dossiers")
        gd.addWidget(self.ed_out, 1, 1)
        b2 = QPushButton("…"); b2.setFixedWidth(28); b2.clicked.connect(self._browse_out)
        gd.addWidget(b2, 1, 2)
        self.btn_scan = QPushButton("🔍  Analyser les RAW")
        self.btn_scan.setObjectName("BtnScan")
        self.btn_scan.clicked.connect(self._scan_raw)
        gd.addWidget(self.btn_scan, 2, 0, 1, 3)
        lv.addWidget(grp_d)

        grp_g = QGroupBox("Regroupement")
        gg = QGridLayout(grp_g)
        gg.addWidget(QLabel("Tolérance (%) :"), 0, 0)
        self.spn_tol = QDoubleSpinBox()
        self.spn_tol.setRange(0.0, 25.0); self.spn_tol.setSingleStep(0.5)
        self.spn_tol.setDecimals(1); self.spn_tol.setValue(3.0)
        self.spn_tol.setToolTip(
            "Deux poses sont mises dans le même groupe si leur écart relatif\n"
            "est sous ce seuil. 3 % fusionne les arrondis (1/1000 vs 1/1024)\n"
            "sans mélanger deux crans de 1/3 de diaph (26 % d'écart).")
        self.spn_tol.valueChanged.connect(self._regroup)
        gg.addWidget(self.spn_tol, 0, 1)
        gg.addWidget(QLabel("Préfixe séquence :"), 1, 0)
        self.ed_prefix = QLineEdit("")
        self.ed_prefix.setPlaceholderText("(vide = nom du dossier)")
        self.ed_prefix.setToolTip("Nom de base des FITS et de la séquence.\n"
                                  "Vide : « 1-160_00001.fit ». Avec « v » : « v1-160_00001.fit ».")
        gg.addWidget(self.ed_prefix, 1, 1)
        lv.addWidget(grp_g)

        grp_c = QGroupBox("Conversion")
        gc = QGridLayout(grp_c)
        self.chk_debayer = QCheckBox("Dématriçage (RAW couleur)")
        self.chk_debayer.setChecked(True)
        self.chk_debayer.setToolTip("Commande Siril « convert -debayer ». "
                                    "Décochez pour garder la matrice de Bayer brute.")
        gc.addWidget(self.chk_debayer, 0, 0, 1, 2)
        gc.addWidget(QLabel("Mise en place :"), 1, 0)
        self.cmb_place = QComboBox()
        self.cmb_place.addItems(["Lien (rapide, 0 octet)", "Copie des RAW"])
        self.cmb_place.setToolTip("Les RAW sont mis à disposition dans le dossier cible\n"
                                  "le temps de la conversion, puis retirés. Le lien évite\n"
                                  "de dupliquer des dizaines de Go.")
        gc.addWidget(self.cmb_place, 1, 1)
        self.chk_32 = QCheckBox("Passer Siril en 32 bits")
        self.chk_32.setToolTip(
            "Réglage global de Siril. Inutile ici (les RAW sont en 14 bits) et\n"
            "double la taille des FITS — mais indispensable AVANT l'empilement\n"
            "pour éviter la postérisation de la couronne faible.")
        gc.addWidget(self.chk_32, 2, 0, 1, 2)
        self.chk_purge = QCheckBox("Vider le dossier cible avant conversion")
        self.chk_purge.setToolTip("Supprime tous les .fit/.fits/.seq du sous-dossier\n"
                                  "cible (pas les RAW). À utiliser en connaissance de cause.")
        gc.addWidget(self.chk_purge, 3, 0, 1, 2)
        lv.addWidget(grp_c)

        hint = QLabel("Un sous-dossier par vitesse est créé dans le dossier de sortie "
                      "(1/160 s → « 1-160 »), avec les FITS dématriçés et la séquence "
                      "prête pour SirilJ Align.")
        hint.setWordWrap(True); hint.setStyleSheet("font-size:8pt; color:#888;")
        lv.addWidget(hint)

        lv.addStretch()
        self.btn_go = QPushButton("▶  Convertir la sélection")
        self.btn_go.setObjectName("BtnGo"); self.btn_go.setEnabled(False)
        self.btn_go.clicked.connect(self._convert)
        lv.addWidget(self.btn_go)
        self.btn_abort = QPushButton("⏹  Arrêter")
        self.btn_abort.setEnabled(False); self.btn_abort.clicked.connect(self._abort_all)
        lv.addWidget(self.btn_abort)
        h.addWidget(left)

        right = QWidget()
        rv = QVBoxLayout(right); rv.setContentsMargins(0, 0, 0, 0); rv.setSpacing(4)
        grp_t = QGroupBox("Plan de rangement")
        gt = QVBoxLayout(grp_t)
        self.table = QTableWidget(0, 4)
        self.table.setHorizontalHeaderLabels(["Vitesse", "Dossier", "Images", "Pose (s)"])
        self.table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        for c in (1, 2, 3):
            self.table.horizontalHeader().setSectionResizeMode(
                c, QHeaderView.ResizeMode.ResizeToContents)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        gt.addWidget(self.table)
        row = QHBoxLayout()
        b_all = QPushButton("Tout cocher"); b_all.clicked.connect(lambda: self._check_all(True))
        b_non = QPushButton("Tout décocher"); b_non.clicked.connect(lambda: self._check_all(False))
        row.addWidget(b_all); row.addWidget(b_non); row.addStretch()
        self.lbl_sum = QLabel("—"); self.lbl_sum.setStyleSheet("color:#999;")
        row.addWidget(self.lbl_sum)
        gt.addLayout(row)
        rv.addWidget(grp_t, stretch=3)
        self.prog = QProgressBar(); rv.addWidget(self.prog)
        self.log = QTextEdit(); self.log.setReadOnly(True); self.log.setFixedHeight(140)
        rv.addWidget(self.log)
        h.addWidget(right, stretch=1)

        self.statusBar().showMessage("Choisissez le dossier des RAW, puis « Analyser ».")

    # =========================================================================
    #  Dossiers
    # =========================================================================
    def _guess_folders(self):
        if self.ed_raw.text().strip():
            return
        cwd = Path(os.getcwd())               # dossier de travail de Siril
        for name in ("raw", "RAW", "Raw"):
            if (cwd / name).is_dir():
                self.ed_raw.setText(str(cwd / name))
                self.ed_out.setText(str(cwd))
                return
        self.ed_raw.setText(str(cwd))
        self.ed_out.setText(str(cwd.parent if cwd.parent != cwd else cwd))

    def _browse_raw(self):
        d = QFileDialog.getExistingDirectory(self, "Dossier des RAW",
                                             self.ed_raw.text() or os.getcwd())
        if d:
            self.ed_raw.setText(d)
            if not self.ed_out.text().strip():
                self.ed_out.setText(str(Path(d).parent))

    def _browse_out(self):
        d = QFileDialog.getExistingDirectory(self, "Dossier de sortie",
                                             self.ed_out.text() or os.getcwd())
        if d:
            self.ed_out.setText(d)

    # =========================================================================
    #  Analyse
    # =========================================================================
    def _scan_raw(self):
        raw_dir = Path(self.ed_raw.text().strip() or ".")
        if not raw_dir.is_dir():
            self._log(f"✗ Dossier introuvable : {raw_dir}")
            return
        files = sorted(str(f) for f in raw_dir.iterdir()
                       if f.is_file() and f.suffix.lower() in RAW_EXT)
        if not files:
            self._log(f"Aucun fichier convertible dans {raw_dir}.")
            return
        self.btn_scan.setEnabled(False); self.btn_go.setEnabled(False)
        self.btn_abort.setEnabled(True); self.prog.setValue(0)
        self._log(f"━━ Analyse de {len(files)} fichier(s) dans {raw_dir}")
        if not shutil.which("exiftool"):
            self._log("(exiftool absent : les RAW récents type CR3 pourraient ne pas "
                      "être lus — installez-le si des poses ressortent « inconnues ».)")
        self._scan = ScanWorker(files)
        self._scan.progress.connect(lambda v, m: (self.prog.setValue(v), self._log(m)))
        self._scan.done.connect(self._on_scan_done)
        self._scan.error.connect(lambda e: self._log(f"✗ {e}"))
        self._scan.finished.connect(self._set_idle)
        self._scan.start()

    def _on_scan_done(self, entries):
        self._entries = entries
        n_ok = sum(1 for e in entries if e.get("sec"))
        self._log(f"✓ Temps de pose lu sur {n_ok}/{len(entries)} fichier(s).")
        self._regroup()

    def _regroup(self):
        if not self._entries:
            return
        self._groups, self._unknown = group_by_speed(self._entries, self.spn_tol.value())
        self._refresh_table()
        if self._unknown:
            self._log(f"⚠ {len(self._unknown)} fichier(s) sans temps de pose lisible — "
                      "non rangés (installez exiftool, ou traitez-les à part).")
        self.btn_go.setEnabled(bool(self._groups))

    def _refresh_table(self):
        self.table.blockSignals(True)
        self.table.setRowCount(len(self._groups))
        for i, g in enumerate(self._groups):
            it0 = QTableWidgetItem(speed_text(g["sec"], g.get("num"), g.get("den")))
            it0.setFlags(it0.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            it0.setCheckState(Qt.CheckState.Checked)
            it1 = QTableWidgetItem(g["label"] + "/")
            it1.setForeground(QColor("#7ec77e"))
            it2 = QTableWidgetItem(str(len(g["files"])))
            it3 = QTableWidgetItem(f"{g['sec']:.6g}")
            for c, itm in enumerate((it0, it1, it2, it3)):
                if c:
                    itm.setFlags(itm.flags() & ~Qt.ItemFlag.ItemIsEditable)
                self.table.setItem(i, c, itm)
        self.table.blockSignals(False)
        self.lbl_sum.setText(f"{len(self._groups)} vitesse(s) · "
                             f"{sum(len(g['files']) for g in self._groups)} image(s)")

    def _check_all(self, state: bool):
        st = Qt.CheckState.Checked if state else Qt.CheckState.Unchecked
        for i in range(self.table.rowCount()):
            it = self.table.item(i, 0)
            if it:
                it.setCheckState(st)

    def _selected_groups(self):
        out = []
        for i, g in enumerate(self._groups):
            it = self.table.item(i, 0)
            if it and it.checkState() == Qt.CheckState.Checked:
                out.append(g)
        return out

    # =========================================================================
    #  Conversion
    # =========================================================================
    def _convert(self):
        groups = self._selected_groups()
        if not groups:
            self._log("Aucun groupe coché.")
            return
        out_root = Path(self.ed_out.text().strip() or ".")
        try:
            out_root.mkdir(parents=True, exist_ok=True)
        except Exception as e:
            self._log(f"✗ Dossier de sortie : {e}")
            return
        params = {
            "groups":     groups,
            "out_root":   str(out_root),
            "prefix":     self.ed_prefix.text().strip(),
            "place_mode": "copie" if self.cmb_place.currentIndex() == 1 else "lien",
            "debayer":    self.chk_debayer.isChecked(),
            "set32":      self.chk_32.isChecked(),
            "purge_all":  self.chk_purge.isChecked(),
        }
        self.btn_go.setEnabled(False); self.btn_scan.setEnabled(False)
        self.btn_abort.setEnabled(True); self.prog.setValue(0)
        self._log(f"━━ Conversion de {len(groups)} vitesse(s) → {out_root}")
        siril_safe_log(self.siril,
                       f"ConvertParVitesse — {len(groups)} vitesse(s) en cours…")
        self._conv = ConvertWorker(self.siril, params)
        self._conv.progress.connect(lambda v, m: (self.prog.setValue(v), self._log(m)))
        self._conv.done.connect(self._on_convert_done)
        self._conv.error.connect(lambda e: self._log(f"✗ {e}"))
        self._conv.finished.connect(self._set_idle)
        self._conv.start()

    def _on_convert_done(self, n_ok, n, aborted):
        self.prog.setValue(100)
        if aborted:
            self._log(f"⏹ Interrompu : {n_ok}/{n} vitesse(s) traitée(s).")
            return
        msg = f"✓ Terminé : {n_ok}/{n} vitesse(s) converties."
        self._log(msg)
        siril_safe_log(self.siril, f"ConvertParVitesse — {msg}")
        if n_ok:
            self._log("→ Dans Siril : Ouvrir une séquence dans le sous-dossier voulu, "
                      "ou lancez SirilJ Align dessus.")

    def _set_idle(self):
        self.btn_scan.setEnabled(True)
        self.btn_abort.setEnabled(False)
        self.btn_go.setEnabled(bool(self._groups))

    def _abort_all(self, wait: bool = False):
        for w in (self._scan, self._conv):
            if w and w.isRunning():
                w.abort()
                if wait:
                    w.wait(5000)
        if not wait:
            self._log("⏹ Arrêt demandé…")

    # =========================================================================
    def _log(self, msg: str):
        esc = html.escape(str(msg))
        weight = "bold" if str(msg).lstrip().startswith(("⚠", "✗")) else "normal"
        self.log.append(f'<div style="font-weight:{weight}">{esc}</div>')
        self.log.ensureCursorVisible()
        self.statusBar().showMessage(str(msg))

    def _restore_settings(self):
        st = QSettings("VeraLux", "ConvertParVitesse")
        if st.contains("geometry"):
            self.restoreGeometry(st.value("geometry"))
        self.ed_raw.setText(st.value("raw_dir", ""))
        self.ed_out.setText(st.value("out_dir", ""))
        self.ed_prefix.setText(st.value("prefix", ""))

    def closeEvent(self, event):
        self._abort_all(wait=True)
        st = QSettings("VeraLux", "ConvertParVitesse")
        st.setValue("geometry", self.saveGeometry())
        st.setValue("raw_dir", self.ed_raw.text())
        st.setValue("out_dir", self.ed_out.text())
        st.setValue("prefix", self.ed_prefix.text())
        super().closeEvent(event)


def main():
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    siril = s.SirilInterface()
    try:
        siril.connect()
    except Exception as e:
        print(f"Avertissement : connexion Siril impossible ({e}).")
    win = ConvertWindow(siril)
    win.show()
    app.exec()


if __name__ == "__main__":
    main()
