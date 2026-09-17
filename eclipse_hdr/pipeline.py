# SPDX-License-Identifier: GPL-3.0-or-later
"""Les étapes du pipeline, câblées sur le projet.

Chaque étape suit le même contrat : elle lit son entrée depuis le projet,
écrit ses sorties dans son dossier, puis marque le manifeste. Elle accepte un
`progress(fait, total, nom)` et un `cancel()` — donc elle est pilotable
indistinctement depuis la CLI ou depuis la GUI, sans rien savoir de l'une ni
de l'autre.

Une étape annulée ne marque pas le manifeste : elle sera reprise entièrement
au prochain lancement, mais les fichiers déjà écrits, eux, sont valides et
seront sautés.
"""

from __future__ import annotations

import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Callable, Optional, Sequence

from . import exposure, imageio, rawdec, stacking
from .project import Project, Timer

__all__ = ["step_scan", "step_decode", "step_stack"]

Progress = Optional[Callable[[int, int, str], None]]
CancelFn = Optional[Callable[[], bool]]


def _noop(*_a, **_k):
    pass


# ── 0. Inventaire ─────────────────────────────────────────────────────────────
def step_scan(project: Project, tol_pct: float = 3.0,
              progress: Progress = None, log=_noop) -> dict:
    """Lit les EXIF de tous les RAW du dossier et les groupe par vitesse.

    Ne décode rien : c'est une étape de quelques secondes, qu'on peut relancer
    librement pour voir le plan avant de s'engager.
    """
    raw_dir = project.raw_dir
    if not raw_dir.is_dir():
        raise FileNotFoundError(f"dossier RAW introuvable : {raw_dir}")

    def prog(i, n, name):
        if progress:
            progress(i, n, name)

    with Timer() as t:
        groups, unknown = exposure.scan_folder(raw_dir, tol_pct, progress=prog)
    project.set_groups(groups, unknown)

    total = sum(len(g["files"]) for g in groups)
    for g in groups:
        log(f"  {g['label']:>8s}  {exposure.speed_text(g['sec'], g.get('num'), g.get('den')):>10s}"
            f"  {len(g['files']):4d} images")
    if unknown:
        log(f"  ⚠ {len(unknown)} fichier(s) sans temps de pose lisible — "
            f"installez exiftool pour les conteneurs récents (CR3…)")

    project.mark("scan", {"tol_pct": tol_pct}, inputs=[raw_dir],
                 extra={"n_groups": len(groups), "n_files": total,
                        "n_unknown": len(unknown)}, elapsed=t.elapsed)
    return {"groups": groups, "unknown": unknown}


# ── 1. Décodage ───────────────────────────────────────────────────────────────
def _decode_one(args):
    """Exécuté dans un processus fils : décode et écrit, ne renvoie que du texte."""
    src, out_tif, out_mask, wb, params, extra, dtype = args
    img, mask, meta = rawdec.decode(src, wb, params)
    meta.update(extra)
    imageio.write(out_tif, img, meta, dtype=dtype)
    imageio.write_mask(out_mask, mask)
    return {"src": Path(src).name, "sat": meta["sat_pixels"],
            "shape": list(img.shape)}


def _wb_key(wb) -> tuple:
    """Signature d'une balance, insensible à un facteur global.

    libraw divise les quatre multiplicateurs par le plus petit avant de les
    appliquer : `2,2,2,2` produit exactement la même image que `1,1,1,1`. Deux
    balances proportionnelles ne sont donc pas un conflit.
    """
    v = [float(x) for x in wb]
    if len(v) == 3:                 # même règle que `_normalize_wb`
        v.append(v[1])
    m = min(x for x in v if x > 0)
    return tuple(round(x / m, 6) for x in v)


def _decoded_labels(project: Project) -> set:
    """Vitesses ayant déjà des images décodées sur le disque."""
    d = project.dir("decode", create=False)
    if not d.is_dir():
        return set()
    return {sub.name for sub in d.iterdir()
            if sub.is_dir() and any(f for f in sub.glob("*.tif")
                                    if not f.name.endswith(".sat.tif"))}


def _setup_reference(project: Project, params: dict, wb_spec=None,
                     force: bool = False, log=_noop) -> tuple[dict, set]:
    """Fixe la balance des blancs et l'orientation, une fois pour la session.

    Ces deux constantes doivent être **identiques pour toutes les images** : la
    balance parce que l'invariant radiométrique en dépend, l'orientation parce
    qu'une image tournée autrement ne s'aligne sur rien.

    Changer la balance alors que des images sont déjà décodées laisserait le
    projet avec deux échelles mélangées — précisément la panne qu'on cherche à
    rendre impossible. C'est donc refusé sans `force` ; avec `force`, toutes les
    vitesses déjà décodées sont refaites, pas seulement celles demandées.

    → (référence, vitesses à refaire pour rester cohérent)
    """
    first = None
    for g in project.groups:
        if g["files"]:
            first = project.abspath(g["files"][0])
            break
    if first is None:
        raise RuntimeError("aucun RAW à décoder — lancez l'inventaire d'abord")

    ref = dict(project.data.get("reference", {}))
    stored = ref.get("wb")

    if wb_spec is None and stored and ref.get("orientation") is not None:
        return ref, set()                     # rien de demandé, rien à changer

    wb, source = rawdec.reference_wb(first, wb_spec)

    stale: set = set()
    if stored and _wb_key(stored) != _wb_key(wb):
        stale = _decoded_labels(project)
        if stale and not force:
            raise ValueError(
                f"la balance des blancs change ({['%.4f' % v for v in stored]} → "
                f"{['%.4f' % v for v in wb]}) alors que "
                f"{len(stale)} vitesse(s) sont déjà décodées.\n"
                f"  Mélanger deux balances dans un même projet casse l'échelle "
                f"radiométrique dont dépend la fusion HDR.\n"
                f"  • pour tout redécoder avec la nouvelle balance : ajoutez --force\n"
                f"  • pour garder l'ancienne : relancez sans --wb")
        if stale:
            log(f"  ⟳ balance changée : {len(stale)} vitesse(s) déjà décodées "
                f"seront refaites ({', '.join(sorted(stale)[:6])}"
                f"{'…' if len(stale) > 6 else ''})")

    if ref.get("orientation") is None:
        orient = rawdec.derive_orientation(first, params.get("demosaic", "AHD"))
        ref["orientation"] = list(orient)
        log(f"  orientation mesurée      : {orient[0] * 90}°"
            f"{' + miroir' if orient[1] else ''} (sur {first.name})")

    ref.update({"wb": wb, "wb_source": source, "measured_on": first.name})
    project.data["reference"] = ref
    project.save()
    log(f"  balance des blancs figée : {['%.4f' % v for v in wb]}  — {source}")
    return ref, stale


def step_decode(project: Project, params: Optional[dict] = None,
                labels: Optional[Sequence[str]] = None, jobs: int = 0,
                force: bool = False, dtype: str = "uint16", wb=None,
                progress: Progress = None, cancel: CancelFn = None,
                log=_noop) -> dict:
    """RAW → TIFF linéaire + masque de saturation, rangés par vitesse.

    `wb` fixe la balance des blancs pour **toute** la session — voir
    `rawdec.reference_wb` pour les formes acceptées. Une fois donnée, elle est
    enregistrée dans le manifeste et réutilisée par les décodages suivants.

    Le format par défaut est **seize bits entiers**, et c'est sans perte : la
    sortie de libraw *est* du seize bits entier, que le float32 se contenterait
    de rembourrer de zéros. Sur une session de quatre cents images, cela fait
    46 Go au lieu de 92. `dtype="float32"` reste disponible pour uniformiser.

    Les étapes suivantes, elles, produisent du float32 : dès l'empilement, la
    moyenne de plusieurs poses crée de vraies valeurs intermédiaires.
    """
    p = dict(rawdec.DecodeParams(**(params or {})))
    ref, stale = _setup_reference(project, p, wb, force, log)
    p["orientation"] = ref["orientation"]
    wb_mul = ref["wb"]

    if stale and labels is not None:
        # La balance a changé : les vitesses déjà décodées doivent suivre,
        # sinon le projet garde deux échelles incompatibles.
        labels = sorted(set(labels) | stale)
    if stale:
        force = True

    groups = [g for g in project.groups
              if labels is None or g["label"] in labels]
    if not groups:
        raise ValueError("aucune vitesse sélectionnée")

    tasks, outputs, skipped = [], [], 0
    for g in groups:
        d = project.dir("decode", g["label"])
        for rel in g["files"]:
            src = project.abspath(rel)
            out_tif = d / (src.stem + ".tif")
            out_mask = d / (src.stem + ".sat.tif")
            outputs.append(out_tif)
            if not force and out_tif.exists() and out_mask.exists():
                skipped += 1
                continue
            extra = {"exptime": g["sec"], "speed": g["label"],
                     "speed_text": exposure.speed_text(g["sec"], g.get("num"),
                                                       g.get("den"))}
            tasks.append((str(src), str(out_tif), str(out_mask), wb_mul, p,
                          extra, dtype))

    if skipped:
        log(f"  {skipped} image(s) déjà décodée(s), ignorée(s)")
    if not tasks:
        log("  rien à faire")
        project.mark("decode", {**p, "dtype": dtype, "wb": wb_mul},
                     inputs=[project.raw_dir], outputs=outputs,
                     extra={"wb": wb_mul, "wb_source": ref.get("wb_source")})
        return {"decoded": 0, "skipped": skipped}

    jobs = jobs or min(os.cpu_count() or 4, 8)
    done, total, sat_total = 0, len(tasks), 0
    with Timer() as t:
        if jobs <= 1:
            for a in tasks:
                if cancel and cancel():
                    log("  ⨯ interrompu")
                    return {"decoded": done, "cancelled": True}
                r = _decode_one(a)
                done += 1
                sat_total += r["sat"]
                if progress:
                    progress(done, total, r["src"])
        else:
            with ProcessPoolExecutor(max_workers=jobs) as ex:
                futures = {ex.submit(_decode_one, a): a for a in tasks}
                try:
                    for f in as_completed(futures):
                        if cancel and cancel():
                            ex.shutdown(wait=False, cancel_futures=True)
                            log("  ⨯ interrompu")
                            return {"decoded": done, "cancelled": True}
                        r = f.result()
                        done += 1
                        sat_total += r["sat"]
                        if progress:
                            progress(done, total, r["src"])
                except BaseException:
                    ex.shutdown(wait=False, cancel_futures=True)
                    raise

    project.mark("decode", {**p, "dtype": dtype, "wb": wb_mul},
                 inputs=[project.raw_dir], outputs=outputs,
                 extra={"wb": wb_mul, "wb_source": ref.get("wb_source"),
                        "orientation": ref["orientation"],
                        "n_decoded": done, "sat_pixels_total": sat_total},
                 elapsed=t.elapsed)
    written = sum(o.stat().st_size for o in outputs if o.exists())
    log(f"  {done} image(s) décodée(s) en {t.elapsed:.0f} s "
        f"({t.elapsed / max(done, 1):.1f} s/image, {jobs} processus), "
        f"{written / 1e9:.1f} Go en {dtype}")
    return {"decoded": done, "skipped": skipped, "seconds": t.elapsed,
            "bytes": written}


# ── 3. Empilement ─────────────────────────────────────────────────────────────
def _stack_source(project: Project) -> str:
    """Empiler les images alignées si elles existent, sinon les décodées."""
    if project.state("align").get("done"):
        return "align"
    return "decode"


def step_stack(project: Project, method: Optional[str] = None,
               low: Optional[float] = None, high: Optional[float] = None,
               missing: str = "nan", labels: Optional[Sequence[str]] = None,
               band_mb: int = 256, force: bool = False,
               progress: Progress = None, cancel: CancelFn = None,
               log=_noop) -> dict:
    """Un master par vitesse.

    Sans méthode explicite, chaque groupe reçoit celle que conseille la
    procédure pour son nombre de poses (médiane ≤ 3, percentile ≤ 6,
    winsorisé au-delà) — le choix est journalisé pour chaque groupe.
    """
    src_step = _stack_source(project)
    groups = [g for g in project.groups
              if labels is None or g["label"] in labels]
    if not groups:
        raise ValueError("aucune vitesse sélectionnée")

    out_dir = project.dir("stack")
    made, outputs, used, absent = [], [], {}, []

    with Timer() as t:
        for g in groups:
            if cancel and cancel():
                log("  ⨯ interrompu")
                return {"masters": len(made), "cancelled": True}

            label = g["label"]
            d = project.dir(src_step, label, create=False)
            frames = sorted((f for f in d.glob("*.tif")
                             if not f.name.endswith(".sat.tif")),
                            key=lambda f: exposure.natural_key(f.name))
            if not frames:
                absent.append(label)      # groupé en une ligne à la fin
                continue

            m, lo, hi = stacking.suggest_method(len(frames))
            if method:
                m = method
                lo = 3.0 if low is None else low
                hi = 3.0 if high is None else high
            else:
                lo = lo if low is None else low
                hi = hi if high is None else high

            out = out_dir / f"master_{label}.tif"
            outputs.append(out)
            if not force and out.exists():
                log(f"  {label:>8s} : master déjà présent, ignoré")
                made.append(out)
                continue

            masks = [f.with_name(f.stem + ".sat.tif") for f in frames]
            masks = masks if all(x.exists() for x in masks) else None

            rej = f" {lo:g}/{hi:g}" if m in ("percentile", "sigma", "winsorized") else ""
            log(f"  {label:>8s} : {len(frames):3d} poses → {m}{rej}")

            def band_progress(y, h, _l=label):
                if progress:
                    progress(y, h, _l)

            r = stacking.stack(
                frames, out, method=m, low=lo, high=hi, missing=missing,
                mask_paths=masks,
                out_mask_path=(out.with_name(f"master_{label}.sat.tif")
                               if masks else None),
                band_bytes=band_mb << 20,
                meta={"exptime": g["sec"], "speed": label,
                      "speed_text": exposure.speed_text(g["sec"], g.get("num"),
                                                        g.get("den")),
                      "source_step": src_step},
                progress=band_progress, cancel=cancel)
            if r is None:
                log("  ⨯ interrompu")
                return {"masters": len(made), "cancelled": True}
            made.append(r)
            used[label] = {"method": m, "low": lo, "high": hi, "n": len(frames)}

    if absent:
        shown = ", ".join(absent[:6]) + ("…" if len(absent) > 6 else "")
        log(f"  {len(absent)} vitesse(s) non décodée(s), ignorée(s) : {shown}")
    project.mark("stack",
                 {"method": method, "low": low, "high": high,
                  "missing": missing, "source": src_step},
                 inputs=[project.dir(src_step, create=False)],
                 outputs=outputs,
                 extra={"per_group": used, "n_masters": len(made),
                        "skipped": absent},
                 elapsed=t.elapsed)
    log(f"  {len(made)} master(s) en {t.elapsed:.0f} s")
    return {"masters": len(made), "per_group": used, "skipped": absent,
            "seconds": t.elapsed}
