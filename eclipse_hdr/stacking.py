# SPDX-License-Identifier: GPL-3.0-or-later
"""Empilement — remplace la commande `stack` de Siril.

Méthodes reprises de Siril/PixInsight, avec les mêmes noms de paramètres pour
que les résultats soient comparables :

| ici | Siril | quand |
|---|---|---|
| `median` | `stack med` | 2 à 3 poses |
| `percentile` | `stack rej p 0.2 0.1` | 4 à 6 poses |
| `winsorized` | `stack rej w 3 3` | 7 poses et plus |
| `sigma` | `stack rej s 3 3` | variante classique |
| `mean` / `sum` / `min` / `max` | idem | cas particuliers |

**Aucune normalisation n'est appliquée, et il n'y a aucun moyen d'en demander
une.** L'équivalent de `-nonorm` est ici le seul comportement possible : toute
normalisation à l'empilement casserait la relation *valeur ∝ pose × luminance*
dont dépend la fusion HDR. Le piège ne peut donc pas se refermer.

**Mémoire bornée.** Les images sont traitées par bandes horizontales : le pic
mémoire dépend du budget demandé, pas du nombre de poses. Empiler quarante
images de 45 Mpx tient dans 300 Mo au lieu de 21 Go.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Callable, Optional, Sequence

import numpy as np
import tifffile

from . import imageio

__all__ = ["METHODS", "combine", "stack"]

METHODS = ("median", "mean", "sum", "min", "max",
           "percentile", "sigma", "winsorized")

_EPS = 1e-9

#: Valeur signalant « pas de donnée ». Un pixel absent est marqué **NaN**, pas
#: zéro : zéro est une mesure légitime (un photosite sous le niveau de noir),
#: et la confondre avec une absence biaise la moyenne vers le haut.
#:
#: Siril, lui, exclut les pixels exactement nuls de ses empilements par
#: moyenne — mais pas de sa médiane. Sur des poses non recalées, cela relève
#: la valeur de quelques pour cent des pixels ; sur des poses recalées, cela
#: protège les bords laissés vides par le décalage. `missing="zero"` reproduit
#: cette convention ; le défaut la refuse, l'alignement d'ici marquant ses
#: bords en NaN.
MISSING = np.float32("nan")


def _nanaware(a: np.ndarray):
    """→ (médiane, écart-type, moyenne) adaptés à la présence de NaN.

    Les variantes `nan*` coûtent trois à cinq fois plus cher : on ne les prend
    que si le bloc contient réellement des valeurs manquantes.
    """
    if np.isnan(a).any():
        return np.nanmedian, np.nanstd, np.nanmean
    return np.median, np.std, np.mean


# ── Combinaison d'un bloc (N, …) ──────────────────────────────────────────────
def _mean_kept(a: np.ndarray, keep: np.ndarray) -> np.ndarray:
    """Moyenne des valeurs conservées ; repli sur la médiane si tout est rejeté."""
    med, _, _ = _nanaware(a)
    cnt = keep.sum(axis=0)
    tot = np.where(keep, a, 0.0).sum(axis=0)      # NaN écarté par `keep`
    out = tot / np.maximum(cnt, 1)
    empty = cnt == 0
    if empty.any():
        fallback = med(a, axis=0)
        out[empty] = np.nan_to_num(fallback, nan=0.0)[empty]
    return out.astype(np.float32)


def _percentile_clip(a, low, high):
    """Écrêtage percentile de Siril : rejet si l'écart à la médiane dépasse la
    fraction demandée **de la médiane**.

    Attention à la nature du critère : il est relatif, donc d'autant plus sévère
    que le signal est faible. Sur le fond de ciel, où la médiane approche zéro,
    il peut tout rejeter — le repli sur la médiane rattrape alors le calcul.
    C'est le comportement de Siril, reproduit tel quel.
    """
    med, _, _ = _nanaware(a)
    m = med(a, axis=0)
    md = np.maximum(np.abs(m), _EPS)
    d = a - m
    with np.errstate(invalid="ignore"):
        keep = ((-d) / md <= low) & (d / md <= high)   # NaN → False, donc écarté
    return _mean_kept(a, keep)


def _sigma_clip(a, low, high, niter=10):
    """Écrêtage sigma classique — centre médian, écart-type ordinaire.

    Reproduit `stack rej s` de Siril, **avec sa faiblesse** : l'écart-type est
    calculé sur des données qui contiennent encore les valeurs aberrantes, donc
    une valeur très éloignée gonfle σ au point de se protéger elle-même. Un
    rayon cosmique à 100 fois le signal n'est pas rejeté.

    La limite est chiffrable : une valeur écartée de *d* parmi *n* poses gonfle
    l'écart-type à √(σ₀² + d²/n), et n'est rejetée à *k* σ que si
    *d*²(1 − k²/n) > k²σ₀². Avec le k = 3 usuel, le facteur s'annule dès que
    **n ≤ 9** : en dessous de dix poses, l'écrêtage sigma ne peut rejeter
    *aucune* valeur isolée, si aberrante soit-elle.

    C'est la raison d'être de la winsorisation, et l'explication de la règle de
    la procédure (percentile sous sept poses). Cette méthode n'est donc jamais
    proposée par `suggest_method` ; elle reste disponible pour reproduire un
    empilement Siril à l'identique.
    """
    keep = ~np.isnan(a)
    with np.errstate(invalid="ignore"):
        for _ in range(niter):
            masked = np.where(keep, a, np.nan)
            m = np.nanmedian(masked, axis=0)
            s = np.nanstd(masked, axis=0)
            new = (a >= m - low * s) & (a <= m + high * s)
            if np.array_equal(new, keep):
                break
            keep = new
    return _mean_kept(a, keep)


def _winsorized_clip(a, low, high, niter=10, tol=5e-4):
    """Écrêtage sigma winsorisé (Huber).

    Sigma est estimé sur les données **winsorisées** à ±1,5 σ — donc insensible
    aux valeurs aberrantes qu'on cherche justement à rejeter — avec le facteur
    correctif 1,134 qui compense le biais introduit par la winsorisation.

    Écart mesuré avec Siril : sur quatre poses, `stack rej w 3 3` de Siril ne
    rejette *rien* — sa sortie est alors identique, bit pour bit, à sa moyenne
    simple. Celle-ci rejette, parce que le sigma winsorisé ne se laisse pas
    gonfler par la valeur aberrante (c'est tout son intérêt) et échappe donc à
    la limite n ≤ k² décrite dans `_sigma_clip`. À vérifier sur vos données :
    la différence ne se manifeste qu'en dessous d'une dizaine de poses.
    """
    med, std, _ = _nanaware(a)
    m = med(a, axis=0)
    s = np.maximum(1.4826 * med(np.abs(a - m), axis=0), _EPS)
    for _ in range(niter):
        w = np.clip(a, m - 1.5 * s, m + 1.5 * s)
        m = med(w, axis=0)
        s_new = np.maximum(1.134 * std(w, axis=0), _EPS)
        done = np.all(np.abs(s_new - s) <= tol * s)
        s = s_new
        if done:
            break
    with np.errstate(invalid="ignore"):
        keep = (a >= m - low * s) & (a <= m + high * s)
    return _mean_kept(a, keep)


def combine(block: np.ndarray, method: str = "winsorized",
            low: float = 3.0, high: float = 3.0,
            missing: str = "nan") -> np.ndarray:
    """Combine un bloc `(N, …)` en `(…)`. N = nombre de poses.

    `missing` dit ce qui compte comme absence de donnée : `"nan"` (défaut) ou
    `"zero"` pour la convention de Siril. Voir `MISSING`.
    """
    a = np.asarray(block, np.float32)
    if missing == "zero":
        a = np.where(a == 0.0, MISSING, a)
    elif missing != "nan":
        raise ValueError(f"convention de donnée manquante inconnue : {missing!r}")
    if a.shape[0] == 1:
        return np.nan_to_num(a[0], nan=0.0)

    nan = np.isnan(a).any()
    if method == "median":
        r = (np.nanmedian if nan else np.median)(a, axis=0)
        return np.nan_to_num(r, nan=0.0).astype(np.float32)
    if method == "mean":
        r = (np.nanmean if nan else np.mean)(a, axis=0)
        return np.nan_to_num(r, nan=0.0).astype(np.float32)
    if method == "sum":
        return (np.nansum if nan else np.sum)(a, axis=0).astype(np.float32)
    if method == "min":
        r = (np.nanmin if nan else np.min)(a, axis=0)
        return np.nan_to_num(r, nan=0.0).astype(np.float32)
    if method == "max":
        r = (np.nanmax if nan else np.max)(a, axis=0)
        return np.nan_to_num(r, nan=0.0).astype(np.float32)
    if method == "percentile":
        return _percentile_clip(a, low, high)
    if method == "sigma":
        return _sigma_clip(a, low, high)
    if method == "winsorized":
        return _winsorized_clip(a, low, high)
    raise ValueError(f"méthode d'empilement inconnue : {method!r} "
                     f"(attendu : {', '.join(METHODS)})")


# ── Empilement de fichiers, par bandes ────────────────────────────────────────
def _open_all(paths):
    """Vues mémoire-mappées + facteur d'échelle de chacune.

    Repli en lecture complète si un TIFF est compressé (la lecture par bandes
    n'est alors pas possible) — l'échelle vaut 1, la lecture l'ayant déjà
    appliquée.
    """
    views, scales, full = [], [], False
    for p in paths:
        try:
            v, s = imageio.open_scaled(p)
            views.append(v)
            scales.append(s)
        except (ValueError, OSError, NotImplementedError):
            views.append(imageio.read(p)[0])
            scales.append(1.0)
            full = True
    return views, scales, full


#: Écart relatif toléré entre temps de pose d'un même empilement — même valeur
#: que la tolérance de regroupement (arrondis APEX du type 1/1000 vs 1/1024).
EXPTIME_TOL = 0.03


def _inherit_meta(paths) -> dict:
    """Métadonnées héritées des images d'entrée, avec contrôle d'homogénéité.

    Empiler des poses de durées différentes détruit l'échelle radiométrique
    sans laisser de trace visible sur le master — l'erreur ne se manifeste que
    plusieurs étapes plus loin, en paliers de luminosité. On la refuse ici,
    au moment où elle est encore diagnosticable.
    """
    metas = []
    for p in paths:
        try:
            metas.append(imageio.read_meta(p))
        except (OSError, ValueError):
            metas.append({})

    exps = [float(m["exptime"]) for m in metas if m.get("exptime")]
    if exps:
        lo, hi = min(exps), max(exps)
        if hi - lo > EXPTIME_TOL * hi:
            raise ValueError(
                f"temps de pose hétérogènes dans l'empilement "
                f"({lo:g} s à {hi:g} s) — chaque vitesse doit être empilée "
                f"séparément, sinon l'échelle radiométrique est perdue")

    out = {}
    for key in ("exptime", "speed", "speed_text", "wb", "orientation",
                "demosaic", "camera"):
        for m in metas:
            if key in m:
                out[key] = m[key]
                break
    return out


def stack(paths: Sequence, out_path, *, method: str = "winsorized",
          low: float = 3.0, high: float = 3.0, missing: str = "nan",
          mask_paths: Optional[Sequence] = None, out_mask_path=None,
          band_bytes: int = 256 << 20, meta: Optional[dict] = None,
          progress: Optional[Callable[[int, int], None]] = None,
          cancel: Optional[Callable[[], bool]] = None) -> Optional[Path]:
    """Empile des TIFF de travail → un master TIFF float32.

    `mask_paths` : masques de saturation des entrées. Le master reçoit la
    **fraction de poses saturées** en chaque pixel — l'information est ainsi
    conservée jusqu'à la fusion HDR au lieu d'être perdue à l'empilement.

    `cancel()` est interrogé à chaque bande : renvoie True pour interrompre.
    Rien n'est alors laissé sur le disque (écriture atomique).

    → chemin du master, ou None si annulé.
    """
    paths = [Path(p) for p in paths]
    if not paths:
        raise ValueError("aucune image à empiler")

    views, scales, full_read = _open_all(paths)
    shapes = {v.shape for v in views}
    if len(shapes) > 1:
        raise ValueError("dimensions hétérogènes : " +
                         ", ".join(f"{Path(p).name}{v.shape}"
                                   for p, v in zip(paths, views)))

    shape = views[0].shape
    h, w = shape[0], shape[1]
    per_px = int(np.prod(shape[1:])) * 4
    n = len(views)
    band = max(1, min(h, band_bytes // max(n * per_px, 1)))

    masks = None
    if mask_paths:
        masks = [imageio.read_mask(m) for m in mask_paths]

    # Le temps de pose est hérité des images d'entrée, pas seulement fourni par
    # l'appelant : c'est la métadonnée dont dépend toute la fusion HDR, et un
    # master qui la perd produit exactement les paliers qu'on cherche à éviter.
    inherited = _inherit_meta(paths)
    m_out = {**inherited, **(meta or {})}
    m_out.update({
        "step": "stack", "linear": True, "method": method,
        "rejection": [low, high] if method in ("percentile", "sigma", "winsorized") else None,
        "stackcnt": n, "missing": missing,
        "sources": [p.name for p in paths],
        "normalization": "none",
    })
    import json
    desc = json.dumps({imageio.META_KEY: imageio.META_VERSION, **m_out},
                      ensure_ascii=False, separators=(",", ":"))

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_path.with_name(out_path.name + ".tmp")
    photometric = "rgb" if len(shape) == 3 and shape[2] == 3 else "minisblack"

    # `ok` couvre les trois sorties possibles — succès, exception, annulation.
    # Sans lui, un `return` d'annulation traverserait le `except` sans rien
    # nettoyer et laisserait justement le fichier tronqué que l'écriture
    # atomique doit rendre impossible.
    ok = False
    try:
        dst = tifffile.memmap(str(tmp), shape=shape, dtype=np.float32,
                              photometric=photometric, description=desc,
                              metadata=None)
        try:
            for y0 in range(0, h, band):
                if cancel and cancel():
                    return None
                y1 = min(y0 + band, h)
                blk = np.empty((n, y1 - y0) + shape[1:], np.float32)
                for i, v in enumerate(views):
                    blk[i] = v[y0:y1]
                    if scales[i] != 1.0:
                        blk[i] *= np.float32(1.0 / scales[i])
                dst[y0:y1] = combine(blk, method, low, high, missing)
                if progress:
                    progress(y1, h)
            dst.flush()
        finally:
            del dst          # démapper avant tout effacement (Windows)
        os.replace(tmp, out_path)
        ok = True
    finally:
        if not ok:
            tmp.unlink(missing_ok=True)
        if not full_read:
            for v in views:
                if isinstance(v, np.memmap):
                    v._mmap.close()          # libère le descripteur tout de suite

    if masks and out_mask_path:
        imageio.write_mask(out_mask_path,
                           np.mean(np.stack(masks, 0), axis=0))
    return out_path


def suggest_method(n: int) -> tuple[str, float, float]:
    """Méthode conseillée selon le nombre de poses — même règle que la procédure.

    L'écrêtage percentile est l'algorithme prévu pour les petits lots ; au-delà
    de six poses, le winsorisé résiste mieux aux valeurs aberrantes (rayons
    cosmiques, pixels chauds, avion ou satellite qui passe).
    """
    if n <= 1:
        return "mean", 0.0, 0.0
    if n <= 3:
        return "median", 0.0, 0.0
    if n <= 6:
        return "percentile", 0.2, 0.1
    return "winsorized", 3.0, 3.0
