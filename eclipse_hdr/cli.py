# SPDX-License-Identifier: GPL-3.0-or-later
"""Ligne de commande : `eclipse-hdr`.

Chaque étape est une sous-commande, exécutable seule et dans n'importe quel
ordre compatible ; `run` les enchaîne. Ctrl-C interrompt proprement : rien
n'est laissé à moitié écrit, et la reprise repart de ce qui existe.
"""

from __future__ import annotations

import argparse
import signal
import sys
from pathlib import Path

import numpy as np

from . import __version__, exposure, imageio, stacking
from .project import Cancel, Project

RAW_SUBDIR = "raw"


# ── Affichage ─────────────────────────────────────────────────────────────────
class Term:
    def __init__(self, quiet=False):
        self.quiet = quiet
        self._bar = False

    def log(self, msg=""):
        if self.quiet:
            return
        if self._bar:
            sys.stdout.write("\n")
            self._bar = False
        print(msg, flush=True)

    def progress(self, done, total, name=""):
        if self.quiet or not total:
            return
        w = 28
        k = int(w * done / total)
        name = (name[:26] + "…") if len(name) > 27 else name
        sys.stdout.write(f"\r  [{'█' * k}{'·' * (w - k)}] {done}/{total}  {name:<28s}")
        sys.stdout.flush()
        self._bar = True
        if done >= total:
            sys.stdout.write("\n")
            self._bar = False

    def done(self):
        if self._bar:
            sys.stdout.write("\n")
            self._bar = False


def _locate(arg: str) -> tuple[Path, Path]:
    """→ (racine du projet, dossier des RAW).

    Le dossier donné est le projet. Ses RAW sont dans `raw/` s'il existe,
    sinon directement dedans — ce qui permet de pointer simplement le dossier
    où sont tous les RAW.
    """
    root = Path(arg).expanduser().resolve()
    sub = root / RAW_SUBDIR
    return (root, sub) if sub.is_dir() else (root, root)


def _speeds(arg):
    return [s.strip() for s in arg.split(",") if s.strip()] if arg else None


# ── Commandes ─────────────────────────────────────────────────────────────────
def cmd_scan(a, term):
    root, raws = _locate(a.dossier)
    if a.raw:
        raws = Path(a.raw).expanduser().resolve()
    if a.project:
        root = Path(a.project).expanduser().resolve()
    if not raws.is_dir():
        term.log(f"✗ dossier RAW introuvable : {raws}")
        return 2

    p = Project.open_or_create(root, raws)
    term.log(f"Projet  : {p.root}")
    term.log(f"RAW     : {raws}")
    from . import pipeline
    pipeline.step_scan(p, a.tol, progress=term.progress, log=term.log)
    term.done()
    n = sum(len(g["files"]) for g in p.groups)
    term.log(f"→ {len(p.groups)} vitesse(s), {n} image(s). Manifeste : {p.root / 'projet.json'}")
    return 0


def cmd_status(a, term):
    p = Project.load(_locate(a.dossier)[0])
    term.log(f"Projet : {p.root}")
    if p.groups:
        term.log("\nVitesses :")
        for g in p.groups:
            term.log(f"  {g['label']:>8s}  "
                     f"{exposure.speed_text(g['sec'], g.get('num'), g.get('den')):>10s}"
                     f"  {len(g['files']):4d} images")
    ref = p.data.get("reference")
    if ref:
        term.log(f"\nBalance des blancs : {['%.4f' % v for v in ref['wb']]}"
                 f"   — {ref.get('wb_source', 'lumière du jour')}")
        term.log(f"Orientation        : {ref['orientation'][0] * 90}°"
                 f"{' + miroir' if ref['orientation'][1] else ''}")
    term.log("\nÉtapes :")
    for s in p.status():
        mark = "✔" if s["done"] else "·"
        det = ""
        if s["done"]:
            det = f"  {s['at'][:16].replace('T', ' ')}"
            if s["n"]:
                det += f"  {s['n']} fichier(s)"
            if s["seconds"]:
                det += f"  {s['seconds']:.0f} s"
        term.log(f"  {mark} {s['label']:<28s}{det}")
    return 0


def cmd_decode(a, term):
    p = Project.load(_locate(a.dossier)[0])
    from . import pipeline
    params = {"demosaic": a.demosaic, "sat_frac": a.sat_frac,
              "sat_dilate": a.sat_dilate, "half_size": a.half}
    cancel = _install_sigint(term)
    r = pipeline.step_decode(p, params, labels=_speeds(a.speeds), jobs=a.jobs,
                             force=a.force, wb=a.wb,
                             dtype="float32" if a.f32 else "uint16",
                             progress=term.progress, cancel=cancel,
                             log=term.log)
    term.done()
    return 1 if r.get("cancelled") else 0


def cmd_stack(a, term):
    p = Project.load(_locate(a.dossier)[0])
    from . import pipeline
    cancel = _install_sigint(term)
    r = pipeline.step_stack(p, method=a.method, low=a.low, high=a.high,
                            missing=a.missing, labels=_speeds(a.speeds),
                            band_mb=a.band,
                            force=a.force, progress=None, cancel=cancel,
                            log=term.log)
    term.done()
    return 1 if r.get("cancelled") else 0


def cmd_run(a, term):
    for name, fn in (("scan", cmd_scan), ("decode", cmd_decode),
                     ("stack", cmd_stack)):
        term.log(f"\n── {name} ─────────────────────────────────────────")
        rc = fn(a, term)
        if rc:
            return rc
    return 0


def cmd_info(a, term):
    import json
    for f in a.fichiers:
        f = Path(f)
        term.log(f"\n{f}")
        try:
            if f.suffix.lower() in exposure.RAW_EXT:
                from . import rawdec
                e = exposure.read_exposure(f)
                if e:
                    term.log(f"  pose : {exposure.speed_text(e[2], e[0], e[1])}"
                             f"  ({e[2]:g} s)")
                for name, wb in rawdec.camera_wb(f).items():
                    lbl = {"daylight": "lumière du jour",
                           "camera": "telle que prise"}[name]
                    term.log(f"  balance {lbl:<16s} : "
                             f"{','.join('%.4f' % v for v in wb)}"
                             f"   → --wb {','.join('%.4f' % v for v in wb[:3])}")
                continue
            if f.suffix.lower() in (".fit", ".fits", ".fts"):
                arr, meta = imageio.from_fits(f)
            else:
                meta = imageio.read_meta(f)
                with __import__("tifffile").TiffFile(str(f)) as tf:
                    pg = tf.pages[0]
                    arr = None
                    term.log(f"  {pg.shape} {pg.dtype}")
            if arr is not None:
                term.log(f"  {arr.shape} {arr.dtype}")
            term.log("  " + json.dumps(meta, ensure_ascii=False, indent=2
                                       ).replace("\n", "\n  "))
        except (OSError, ValueError) as e:
            term.log(f"  ✗ {e}")
    return 0


def cmd_export(a, term):
    src = Path(a.fichier)
    arr, meta = (imageio.from_fits(src) if src.suffix.lower() in
                 (".fit", ".fits", ".fts") else imageio.read(src))
    out = Path(a.out) if a.out else src.with_name(src.stem + "_16b.tif")
    imageio.export16(out, arr, transform=a.transform, asinh_k=a.asinh_k,
                     scale=a.scale, meta=meta)
    term.log(f"→ {out}  ({a.transform}, 16 bits)")
    return 0


def cmd_diff(a, term):
    """Compare deux images — l'outil de validation contre Siril."""
    def load(p):
        p = Path(p)
        return (imageio.from_fits(p)[0] if p.suffix.lower() in
                (".fit", ".fits", ".fts") else imageio.read(p)[0])

    x, y = load(a.a), load(a.b)
    if x.shape != y.shape:
        term.log(f"✗ dimensions différentes : {x.shape} vs {y.shape}")
        return 2
    if a.normalise:
        sx, sy = float(np.median(x)), float(np.median(y))
        if sy > 0:
            y = y * (sx / sy)
            term.log(f"  (b remise à l'échelle ×{sx / sy:.6g} sur la médiane)")
    d = x.astype(np.float64) - y.astype(np.float64)
    scale = float(np.percentile(np.abs(x), 99.9)) or 1.0
    term.log(f"  a : min {x.min():.6g}  méd {np.median(x):.6g}  max {x.max():.6g}")
    term.log(f"  b : min {y.min():.6g}  méd {np.median(y):.6g}  max {y.max():.6g}")
    term.log(f"  écart  RMS {np.sqrt((d ** 2).mean()):.6g}  "
             f"max |Δ| {np.abs(d).max():.6g}")
    term.log(f"  relatif au 99,9ᵉ centile de a ({scale:.6g}) : "
             f"RMS {np.sqrt((d ** 2).mean()) / scale * 100:.4f} %  "
             f"max {np.abs(d).max() / scale * 100:.4f} %")
    return 0


# ── Interruption ──────────────────────────────────────────────────────────────
def _install_sigint(term) -> Cancel:
    """Premier Ctrl-C : arrêt propre après le fichier en cours. Second : brutal."""
    c = Cancel()

    def handler(_sig, _frm):
        if c:
            raise KeyboardInterrupt
        c.set()
        term.done()
        term.log("  ⏸  arrêt demandé — fin du fichier en cours "
                 "(Ctrl-C à nouveau pour forcer)")

    signal.signal(signal.SIGINT, handler)
    return c


# ── Analyseur ─────────────────────────────────────────────────────────────────
def build_parser():
    ap = argparse.ArgumentParser(
        prog="eclipse-hdr",
        description="Traitement HDR d'éclipse solaire, du RAW à l'image finale. "
                    "Chaque étape produit un TIFF 32 bits linéaire et peut être "
                    "lancée, interrompue et reprise séparément.")
    ap.add_argument("--version", action="version", version=f"eclipse-hdr {__version__}")
    ap.add_argument("-q", "--quiet", action="store_true")
    sub = ap.add_subparsers(dest="cmd", required=True)

    def add_dir(sp):
        sp.add_argument("dossier", help="dossier du projet (contenant les RAW, "
                                        "ou un sous-dossier raw/)")

    s = sub.add_parser("scan", help="inventorier les RAW et les grouper par vitesse")
    add_dir(s)
    s.add_argument("--raw", help="dossier des RAW s'il est ailleurs")
    s.add_argument("--project", help="écrire le projet ailleurs")
    s.add_argument("--tol", type=float, default=3.0,
                   help="tolérance de regroupement, en %% (défaut : 3)")
    s.set_defaults(func=cmd_scan)

    s = sub.add_parser("status", help="état du projet")
    add_dir(s)
    s.set_defaults(func=cmd_status)

    s = sub.add_parser("decode", help="RAW → TIFF 32 bits linéaire + masque de saturation")
    add_dir(s)
    s.add_argument("--speeds", help="vitesses à traiter, séparées par des virgules")
    s.add_argument("--jobs", type=int, default=0, help="processus (défaut : cœurs, max 8)")
    s.add_argument("--demosaic", default="AHD", choices=sorted(__import__(
        "eclipse_hdr.rawdec", fromlist=["DEMOSAIC"]).DEMOSAIC))
    s.add_argument("--sat-frac", dest="sat_frac", type=float, default=0.98,
                   help="seuil de saturation, en fraction de la plage utile")
    s.add_argument("--sat-dilate", dest="sat_dilate", type=int, default=3,
                   help="dilatation du masque, en px (portée du dématriçage)")
    s.add_argument("--wb", metavar="BALANCE", help="balance des blancs, la même pour toutes les images : « daylight » (défaut, constante du boîtier), « camera » (telle que prise), « neutral », des valeurs « R,V,B » ou « R,V,B,V2 », ou le chemin d'un RAW dont reprendre la balance. Seuls les rapports comptent : 2,2,2,2 est identique à 1,1,1,1")
    s.add_argument("--half", action="store_true", help="quart de résolution (aperçu)")
    s.add_argument("--f32", action="store_true",
                   help="écrire en float32 plutôt qu'en 16 bits entiers "
                        "(deux fois plus de disque, aucune information de plus : "
                        "la sortie de libraw est déjà du 16 bits)")
    s.add_argument("--force", action="store_true", help="refaire même si présent")
    s.set_defaults(func=cmd_decode)

    s = sub.add_parser("stack", help="empiler chaque vitesse en un master")
    add_dir(s)
    s.add_argument("--speeds")
    s.add_argument("--method", choices=stacking.METHODS,
                   help="défaut : choisie par groupe selon le nombre de poses")
    s.add_argument("--low", type=float, help="seuil de rejet bas")
    s.add_argument("--high", type=float, help="seuil de rejet haut")
    s.add_argument("--band", type=int, default=256,
                   help="budget mémoire par bande, en Mo (défaut : 256)")
    s.add_argument("--missing", default="nan", choices=("nan", "zero"),
                   help="ce qui compte comme absence de donnée. « zero » "
                        "reproduit la convention de Siril, qui exclut les "
                        "pixels exactement nuls de ses moyennes")
    s.add_argument("--force", action="store_true")
    s.set_defaults(func=cmd_stack)

    s = sub.add_parser("run", help="enchaîner les étapes disponibles")
    add_dir(s)
    s.add_argument("--raw")
    s.add_argument("--project")
    s.add_argument("--tol", type=float, default=3.0)
    s.add_argument("--speeds")
    s.add_argument("--jobs", type=int, default=0)
    s.add_argument("--demosaic", default="AHD")
    s.add_argument("--sat-frac", dest="sat_frac", type=float, default=0.98)
    s.add_argument("--sat-dilate", dest="sat_dilate", type=int, default=3)
    s.add_argument("--wb", metavar="BALANCE")
    s.add_argument("--half", action="store_true")
    s.add_argument("--f32", action="store_true")
    s.add_argument("--method", choices=stacking.METHODS)
    s.add_argument("--low", type=float)
    s.add_argument("--high", type=float)
    s.add_argument("--band", type=int, default=256)
    s.add_argument("--missing", default="nan", choices=("nan", "zero"))
    s.add_argument("--force", action="store_true")
    s.set_defaults(func=cmd_run)

    s = sub.add_parser("info", help="afficher les métadonnées d'une image")
    s.add_argument("fichiers", nargs="+")
    s.set_defaults(func=cmd_info)

    s = sub.add_parser("export", help="exporter en TIFF 16 bits (Photoshop, DxO…)")
    s.add_argument("fichier")
    s.add_argument("-o", "--out")
    s.add_argument("--transform", default="autostretch",
                   choices=("autostretch", "asinh", "linear"))
    s.add_argument("--asinh-k", dest="asinh_k", type=float, default=100.0)
    s.add_argument("--scale", type=float, help="facteur pour --transform linear")
    s.set_defaults(func=cmd_export)

    s = sub.add_parser("diff", help="comparer deux images (validation)")
    s.add_argument("a")
    s.add_argument("b")
    s.add_argument("--normalise", action="store_true",
                   help="remettre b à l'échelle de a sur la médiane avant de comparer")
    s.set_defaults(func=cmd_diff)

    return ap


def main(argv=None) -> int:
    a = build_parser().parse_args(argv)
    term = Term(a.quiet)
    try:
        return a.func(a, term)
    except KeyboardInterrupt:
        term.done()
        term.log("⨯ interrompu")
        return 130
    except (FileNotFoundError, ValueError, KeyError, RuntimeError) as e:
        term.done()
        term.log(f"✗ {e}")
        return 2


if __name__ == "__main__":
    sys.exit(main())
