# SPDX-License-Identifier: GPL-3.0-or-later
"""Projet et manifeste — l'état du traitement, sur le disque.

Un traitement est un dossier auto-descriptif : les RAW, un manifeste JSON, et
un sous-dossier par étape. Le manifeste porte les groupes de vitesse, les
paramètres de chaque étape et son état. Il rend le traitement **reprenable**
(on sait ce qui est fait) et **reproductible** (on sait avec quels réglages).

Reprise incrémentale
--------------------
Chaque étape enregistre une **empreinte** : `sha256(nom + paramètres +
empreintes des entrées)`. Relancer le pipeline ne recalcule que les étapes
dont l'empreinte a changé. Retoucher un réglage de la couronne ne réempile
donc pas les poses. C'est un système de dépendances minimal, volontairement :
une chaîne linéaire, pas un graphe.

L'empreinte d'un fichier est `(taille, date de modification)`, pas son
contenu : hacher des masters de 500 Mo coûterait plus cher que de les
recalculer.
"""

from __future__ import annotations

import hashlib
import json
import os
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Optional

__all__ = ["Project", "Cancel", "STEPS", "STEP_DIRS"]

MANIFEST = "projet.json"
SCHEMA = 1

# Ordre du pipeline. `None` = étape sans sortie image (inventaire).
STEPS = ["scan", "decode", "align", "stack", "coalign",
         "hdr", "corona", "composite", "export"]

STEP_DIRS = {
    "decode": "01_decode",
    "align": "02_align",
    "stack": "03_stack",
    "coalign": "04_coalign",
    "hdr": "05_hdr",
    "corona": "06_corona",
    "composite": "07_composite",
    "export": "export",
}

STEP_LABELS = {
    "scan": "Inventaire des RAW",
    "decode": "Décodage RAW",
    "align": "Alignement sur la Lune",
    "stack": "Empilement par vitesse",
    "coalign": "Co-alignement des masters",
    "hdr": "Fusion HDR",
    "corona": "Révélation de la couronne",
    "composite": "Composition par couches",
    "export": "Export 16 bits",
}


class Cancel:
    """Jeton d'annulation coopératif, interrogé par les étapes longues."""

    def __init__(self):
        self._flag = False

    def set(self):
        self._flag = True

    def clear(self):
        self._flag = False

    def __bool__(self):
        return self._flag

    def __call__(self):
        return self._flag


def _now() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


def file_stamp(path) -> str:
    """Empreinte bon marché d'un fichier : taille + date de modification."""
    try:
        st = os.stat(path)
    except OSError:
        return "absent"
    return f"{st.st_size}:{st.st_mtime_ns}"


def _canon(obj: Any) -> str:
    return json.dumps(obj, sort_keys=True, ensure_ascii=False,
                      separators=(",", ":"), default=str)


class Project:
    """Un dossier de traitement."""

    def __init__(self, root, data: Optional[dict] = None):
        self.root = Path(root).resolve()
        self.data = data if data is not None else self._blank()

    # ── Création / chargement ────────────────────────────────────────────
    @staticmethod
    def _blank() -> dict:
        return {"eclipse_hdr_project": SCHEMA, "created": _now(),
                "raw_dir": "raw", "groups": [], "unknown": [], "steps": {}}

    @classmethod
    def create(cls, root, raw_dir=None) -> "Project":
        root = Path(root).resolve()
        root.mkdir(parents=True, exist_ok=True)
        p = cls(root)
        if raw_dir:
            rd = Path(raw_dir).resolve()
            try:
                p.data["raw_dir"] = str(rd.relative_to(root))
            except ValueError:
                p.data["raw_dir"] = str(rd)      # RAW hors du projet : absolu
        p.save()
        return p

    @classmethod
    def load(cls, root) -> "Project":
        root = Path(root).resolve()
        f = root / MANIFEST
        if not f.exists():
            raise FileNotFoundError(
                f"pas de projet dans {root} (aucun {MANIFEST}). "
                f"Lancez d'abord « eclipse-hdr scan ».")
        with open(f, encoding="utf-8") as fh:
            data = json.load(fh)
        ver = data.get("eclipse_hdr_project")
        if ver != SCHEMA:
            raise ValueError(f"manifeste en version {ver}, attendu {SCHEMA}")
        return cls(root, data)

    @classmethod
    def open_or_create(cls, root, raw_dir=None) -> "Project":
        try:
            return cls.load(root)
        except FileNotFoundError:
            return cls.create(root, raw_dir)

    def save(self) -> Path:
        """Écriture atomique — un manifeste n'est jamais à moitié écrit."""
        f = self.root / MANIFEST
        tmp = f.with_name(f.name + ".tmp")
        self.data["saved"] = _now()
        try:
            with open(tmp, "w", encoding="utf-8") as fh:
                json.dump(self.data, fh, ensure_ascii=False, indent=1)
            os.replace(tmp, f)
        except BaseException:
            tmp.unlink(missing_ok=True)
            raise
        return f

    # ── Chemins ──────────────────────────────────────────────────────────
    @property
    def raw_dir(self) -> Path:
        rd = Path(self.data.get("raw_dir", "raw"))
        return rd if rd.is_absolute() else self.root / rd

    def dir(self, step: str, group: Optional[str] = None,
            create: bool = True) -> Path:
        if step not in STEP_DIRS:
            raise KeyError(f"étape sans dossier : {step!r}")
        d = self.root / STEP_DIRS[step]
        if group:
            d = d / group
        if create:
            d.mkdir(parents=True, exist_ok=True)
        return d

    # ── Groupes de vitesse ───────────────────────────────────────────────
    @property
    def groups(self) -> list[dict]:
        return self.data.get("groups", [])

    def group(self, label: str) -> dict:
        for g in self.groups:
            if g["label"] == label:
                return g
        raise KeyError(f"vitesse inconnue : {label!r} "
                       f"(connues : {', '.join(g['label'] for g in self.groups)})")

    def set_groups(self, groups: Iterable[dict], unknown: Iterable[dict] = ()):
        """Enregistre l'inventaire, en chemins relatifs au projet si possible."""
        def rel(p):
            try:
                return str(Path(p).resolve().relative_to(self.root))
            except ValueError:
                return str(p)

        self.data["groups"] = [{
            "label": g["label"], "sec": g["sec"],
            "num": g.get("num"), "den": g.get("den"),
            "files": [rel(e["path"]) for e in g["files"]],
        } for g in groups]
        self.data["unknown"] = [rel(e["path"]) for e in unknown]

    def abspath(self, p) -> Path:
        p = Path(p)
        return p if p.is_absolute() else self.root / p

    # ── État des étapes ──────────────────────────────────────────────────
    def fingerprint(self, step: str, params: dict,
                    inputs: Iterable = ()) -> str:
        h = hashlib.sha256()
        h.update(step.encode())
        h.update(_canon(params).encode())
        for i in sorted(str(x) for x in inputs):
            h.update(i.encode())
            h.update(file_stamp(i).encode())
        return h.hexdigest()[:16]

    def state(self, step: str) -> dict:
        return self.data.setdefault("steps", {}).get(step, {})

    def is_current(self, step: str, params: dict, inputs: Iterable = (),
                   outputs: Iterable = ()) -> bool:
        """L'étape est-elle à jour ? Faux si les paramètres ou les entrées ont
        changé, ou si une sortie attendue manque."""
        st = self.state(step)
        if not st.get("done"):
            return False
        if st.get("fp") != self.fingerprint(step, params, inputs):
            return False
        return all(Path(o).exists() for o in outputs)

    def mark(self, step: str, params: dict, inputs: Iterable = (),
             outputs: Iterable = (), extra: Optional[dict] = None,
             elapsed: Optional[float] = None):
        """Marque une étape terminée et enregistre le manifeste."""
        def rel(p):
            try:
                return str(Path(p).resolve().relative_to(self.root))
            except ValueError:
                return str(p)

        outs = [rel(o) for o in outputs]
        self.data.setdefault("steps", {})[step] = {
            "done": True, "at": _now(),
            "fp": self.fingerprint(step, params, inputs),
            "params": params, "outputs": outs, "n_outputs": len(outs),
            **({"seconds": round(elapsed, 1)} if elapsed is not None else {}),
            **(extra or {}),
        }
        self.save()

    def invalidate(self, step: str, cascade: bool = True):
        """Marque une étape (et les suivantes) comme à refaire."""
        steps = self.data.setdefault("steps", {})
        idx = STEPS.index(step)
        for s in STEPS[idx:] if cascade else [step]:
            steps.pop(s, None)
        self.save()

    def status(self) -> list[dict]:
        """Résumé pour l'affichage : une ligne par étape."""
        out = []
        for s in STEPS:
            st = self.state(s)
            out.append({"step": s, "label": STEP_LABELS[s],
                        "done": bool(st.get("done")), "at": st.get("at"),
                        "n": st.get("n_outputs"), "seconds": st.get("seconds")})
        return out


class Timer:
    """Chronomètre de bloc, pour renseigner `elapsed` dans le manifeste."""

    def __enter__(self):
        self.t0 = time.monotonic()
        return self

    def __exit__(self, *exc):
        self.elapsed = time.monotonic() - self.t0
        return False
