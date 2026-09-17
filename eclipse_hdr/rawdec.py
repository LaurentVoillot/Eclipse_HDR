# SPDX-License-Identifier: GPL-3.0-or-later
"""Décodage RAW linéaire, avec masque de saturation lu dans le domaine capteur.

Remplace `convert -debayer` de Siril par libraw (via rawpy), avec trois
garanties que la conversion Siril ne donne pas :

**1. L'invariant radiométrique.** `gamma=(1,1)`, `no_auto_bright=True` et des
multiplicateurs de balance des blancs **figés, identiques pour toutes les
images**. Sans cela chaque pose reçoit son propre facteur d'échelle et la
relation *valeur ∝ pose × luminance* — dont dépend toute la fusion HDR — est
détruite. C'est la cause racine des paliers de luminosité.

**2. La saturation exacte.** Après balance des blancs et matrice couleur, un
photosite saturé n'atterrit plus à 65535 : chaque canal est multiplié
différemment, puis mélangé par la matrice. Détecter la saturation sur la
sortie ne peut donc être qu'approximatif. On la lit ici **avant dématriçage**,
par comparaison au niveau de saturation du capteur, puis on la dilate du rayon
du dématriçage — un photosite saturé contamine ses voisins interpolés.

**3. Pas de reconstruction des hautes lumières.** `highlight_mode=Clip` : on
veut que la saturation reste franche et détectable, pas que libraw invente des
valeurs plausibles là où l'information est perdue.

La rotation appliquée par l'appareil est **mesurée**, pas devinée : voir
`derive_orientation`.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
import rawpy

__all__ = ["DEMOSAIC", "reference_wb", "camera_wb", "derive_orientation",
           "decode", "DecodeParams"]

DEMOSAIC = {
    "AHD": rawpy.DemosaicAlgorithm.AHD,
    "VNG": rawpy.DemosaicAlgorithm.VNG,
    "PPG": rawpy.DemosaicAlgorithm.PPG,
    "DHT": rawpy.DemosaicAlgorithm.DHT,
    "AAHD": rawpy.DemosaicAlgorithm.AAHD,
    "LINEAR": rawpy.DemosaicAlgorithm.LINEAR,
}


class DecodeParams(dict):
    """Paramètres de décodage, sérialisables dans le manifeste."""

    DEFAULTS = {
        "demosaic": "AHD",
        "wb": None,            # None → balance lumière du jour de l'appareil
        "sat_frac": 0.98,      # seuil de saturation, en fraction de la plage utile
        "sat_dilate": 3,       # px — portée du dématriçage
        "orientation": None,   # None → mesurée sur le premier fichier
        "half_size": False,    # aperçu rapide (quart de résolution)
    }

    def __init__(self, **kw):
        super().__init__({**self.DEFAULTS, **kw})


# ── Balance des blancs de référence ───────────────────────────────────────────
def _normalize_wb(wb, where: str) -> list[float]:
    """Complète à quatre multiplicateurs et vérifie qu'ils sont exploitables."""
    wb = [float(v) for v in wb]
    if len(wb) == 3:
        wb.append(wb[1])            # capteur RGBG : le second vert suit le premier
    if len(wb) < 4:
        wb += [0.0] * (4 - len(wb))
    wb = wb[:4]
    if wb[3] == 0.0:
        wb[3] = wb[1]
    if any(v < 0 for v in wb) or not any(v > 0 for v in wb[:3]):
        raise ValueError(f"balance des blancs inexploitable ({where}) : {wb}")
    return wb


def camera_wb(path) -> dict[str, list[float]]:
    """Les deux balances inscrites dans un RAW, pour choisir en connaissance."""
    with rawpy.imread(str(path)) as raw:
        out = {}
        for name, val in (("daylight", raw.daylight_whitebalance),
                          ("camera", raw.camera_whitebalance)):
            v = list(map(float, val))
            if any(x > 0 for x in v):
                out[name] = _normalize_wb(v, f"{name} de {Path(path).name}")
    return out


def reference_wb(path, spec=None) -> tuple[list[float], str]:
    """Multiplicateurs à appliquer **à toutes** les images de la session.

    `spec` accepte :

    | valeur | balance retenue |
    |---|---|
    | `None` ou `"daylight"` | « lumière du jour » de l'appareil *(défaut)* |
    | `"camera"` | telle que prise, lue sur `path` |
    | `"neutral"` | aucune : multiplicateurs égaux |
    | `"1.97,0.94,1.14"` | valeurs explicites (3 ou 4) |
    | un chemin de RAW | balance telle que prise de **ce** fichier |

    Le défaut est la balance « lumière du jour » parce que c'est une **constante
    du boîtier** : elle vaut la même chose pour tous les fichiers, même si la
    prise de vue était en balance automatique. Mesuré sur vingt NEF, la balance
    *telle que prise* prenait vingt valeurs différentes — la retenir donnerait à
    chaque image son propre facteur d'échelle et détruirait l'invariant
    radiométrique dont dépend la fusion HDR.

    **Seuls les rapports comptent** : libraw divise les quatre multiplicateurs
    par le plus petit avant de les appliquer. `2,2,2,2` est donc identique à
    `1,1,1,1` — ce réglage change la teinte, jamais la luminosité.

    → (multiplicateurs, description de la provenance)
    """
    if spec is None or spec == "" or spec == "daylight":
        found = camera_wb(path)
        if "daylight" in found:
            return found["daylight"], f"lumière du jour ({Path(path).name})"
        if "camera" in found:
            return found["camera"], (f"telle que prise ({Path(path).name}) — "
                                     f"pas de balance lumière du jour dans ce RAW")
        raise ValueError(f"aucune balance des blancs exploitable dans {path}")

    if spec == "neutral" or spec == "none":
        return [1.0, 1.0, 1.0, 1.0], "neutre (aucune balance)"

    if spec == "camera":
        found = camera_wb(path)
        if "camera" not in found:
            raise ValueError(f"pas de balance « telle que prise » dans {path}")
        return found["camera"], f"telle que prise ({Path(path).name})"

    if "," in str(spec):
        try:
            vals = [float(v) for v in str(spec).replace(" ", "").split(",") if v]
        except ValueError:
            raise ValueError(
                f"balance des blancs illisible : {spec!r} — attendu trois ou "
                f"quatre nombres séparés par des virgules, p. ex. « 1.97,0.94,1.14 »"
            ) from None
        if len(vals) not in (3, 4):
            raise ValueError(f"balance des blancs : {len(vals)} valeur(s) données, "
                             f"il en faut trois (RVB) ou quatre (RVBV)")
        return _normalize_wb(vals, "valeurs explicites"), "valeurs explicites"

    p = Path(str(spec)).expanduser()
    if p.is_file():
        found = camera_wb(p)
        if "camera" not in found:
            raise ValueError(f"pas de balance « telle que prise » dans {p}")
        return found["camera"], f"telle que prise ({p.name})"

    raise ValueError(
        f"balance des blancs inconnue : {spec!r} — attendu « daylight », "
        f"« camera », « neutral », des valeurs « R,V,B », ou un chemin de RAW")


# ── Orientation ───────────────────────────────────────────────────────────────
_DIHEDRAL = [("rot", 0, False), ("rot", 1, False), ("rot", 2, False),
             ("rot", 3, False), ("rot", 0, True), ("rot", 1, True),
             ("rot", 2, True), ("rot", 3, True)]


def _apply_orientation(a: np.ndarray, orient) -> np.ndarray:
    """Applique (k quarts de tour, miroir) — identiquement à l'image et au masque."""
    if not orient:
        return a
    k, mirror = orient
    out = np.rot90(a, k) if k else a
    if mirror:
        out = out[:, ::-1]
    return np.ascontiguousarray(out)


def derive_orientation(path, demosaic: str = "AHD") -> tuple[int, bool]:
    """Mesure la rotation que l'appareil demande, → (quarts de tour, miroir).

    On décode deux fois le **même** fichier en quart de résolution — une fois
    sans rotation, une fois avec celle de l'appareil — puis on cherche laquelle
    des huit transformations du carré fait coïncider les deux. C'est plus sûr
    que de réinterpréter le code `flip` de libraw, dont la sémantique dépend
    d'un ordre transposition/miroir facile à inverser par erreur.

    Une seule mesure par session : l'orientation est une constante du boîtier.
    """
    kw = dict(gamma=(1, 1), no_auto_bright=True, output_bps=16, half_size=True,
              use_camera_wb=True, use_auto_wb=False,
              demosaic_algorithm=DEMOSAIC[demosaic],
              highlight_mode=rawpy.HighlightMode.Clip)
    with rawpy.imread(str(path)) as raw:
        plain = raw.postprocess(user_flip=0, **kw)
    with rawpy.imread(str(path)) as raw:
        cam = raw.postprocess(**kw)                 # rotation de l'appareil

    p = plain[::4, ::4].astype(np.float32)
    c = cam[::4, ::4].astype(np.float32)
    best, best_err = (0, False), np.inf
    for _, k, mirror in _DIHEDRAL:
        t = _apply_orientation(p, (k, mirror))
        if t.shape != c.shape:
            continue
        err = float(np.mean(np.abs(t - c)))
        if err < best_err:
            best, best_err = (k, mirror), err

    scale = float(np.mean(np.abs(c))) or 1.0
    if best_err > 0.02 * scale:      # aucune transformation ne colle
        return (0, False)
    return best


# ── Masque de saturation ──────────────────────────────────────────────────────
def _saturation_mask(raw, sat_frac: float, dilate: int) -> np.ndarray:
    """Masque booléen (H,W) dans la géométrie capteur, avant dématriçage."""
    visible = raw.raw_image_visible
    colors = raw.raw_colors_visible

    black = np.asarray(raw.black_level_per_channel, dtype=np.float32)
    cw = getattr(raw, "camera_white_level_per_channel", None)
    if cw and any(v for v in cw):
        white = np.asarray(cw, dtype=np.float32)
    else:
        white = np.full(4, float(raw.white_level), dtype=np.float32)
    if black.size < 4:
        black = np.resize(black, 4)
    if white.size < 4:
        white = np.resize(white, 4)

    thr = black + (white - black) * float(sat_frac)
    sat = visible >= thr[colors]

    if dilate and dilate > 0:
        import cv2
        k = 2 * int(dilate) + 1
        sat = cv2.dilate(sat.astype(np.uint8), np.ones((k, k), np.uint8)) > 0
    return sat


def _match_shape(mask: np.ndarray, hw: tuple[int, int]) -> np.ndarray:
    """Recadre/complète le masque sur la géométrie de sortie de libraw.

    Les deux coïncident sur la quasi-totalité des boîtiers ; quand elles
    diffèrent (marges de recadrage), on recadre au centre plutôt que de
    renoncer au masque.
    """
    h, w = hw
    mh, mw = mask.shape
    if (mh, mw) == (h, w):
        return mask
    out = np.zeros((h, w), bool)
    ch, cw_ = min(mh, h), min(mw, w)
    sy, sx = (mh - ch) // 2, (mw - cw_) // 2
    dy, dx = (h - ch) // 2, (w - cw_) // 2
    out[dy:dy + ch, dx:dx + cw_] = mask[sy:sy + ch, sx:sx + cw_]
    return out


# ── Décodage ──────────────────────────────────────────────────────────────────
def decode(path, wb: list[float], params: Optional[dict] = None):
    """RAW → (image (H,W,3) float32 linéaire dans [0,1], masque bool, méta).

    `wb` doit être **le même pour toute la session** (voir `reference_wb`).
    """
    p = DecodeParams(**(params or {}))
    path = Path(path)

    with rawpy.imread(str(path)) as raw:
        mask = _saturation_mask(raw, p["sat_frac"], p["sat_dilate"])
        if p["half_size"]:
            mask = mask[::2, ::2]
        rgb = raw.postprocess(
            user_flip=0,                       # rotation appliquée par nos soins
            gamma=(1, 1),                      # linéaire
            no_auto_bright=True,               # aucune mise à l'échelle par image
            use_camera_wb=False, use_auto_wb=False, user_wb=list(wb),
            output_bps=16,
            output_color=rawpy.ColorSpace.sRGB,     # matrice 3×3 : reste linéaire
            demosaic_algorithm=DEMOSAIC[p["demosaic"]],
            highlight_mode=rawpy.HighlightMode.Clip,
            median_filter_passes=0,
            fbdd_noise_reduction=rawpy.FBDDNoiseReductionMode.Off,
            half_size=bool(p["half_size"]),
            four_color_rgb=False,
        )
        black = list(map(float, raw.black_level_per_channel))
        white = float(raw.white_level)

    mask = _match_shape(mask, rgb.shape[:2])

    orient = p["orientation"]
    if orient:
        orient = (int(orient[0]), bool(orient[1]))
        rgb = _apply_orientation(rgb, orient)
        mask = _apply_orientation(mask, orient)

    img = rgb.astype(np.float32)
    img *= np.float32(1.0 / 65535.0)

    meta = {
        "step": "decode", "linear": True, "source": path.name,
        "wb": [float(v) for v in wb], "demosaic": p["demosaic"],
        "sat_frac": p["sat_frac"], "sat_dilate": p["sat_dilate"],
        "orientation": list(orient) if orient else [0, False],
        "black_level": black, "white_level": white,
        "sat_pixels": int(mask.sum()),
    }
    return img, mask, meta
