# Eclipse_HDR

Outils Python pour [Siril](https://siril.org/) dédiés au **traitement HDR d'éclipse solaire** — et réutilisables pour le lunaire, le planétaire et le ciel profond.

*Python tools for Siril focused on solar-eclipse HDR processing — also usable for lunar, planetary and deep-sky imaging. The UI is in French; bilingual guides are included.*

Les scripts s'exécutent depuis **Siril → Scripts**, ouvrent une interface PyQt6 (thème sombre) et installent leurs dépendances au premier lancement via `sirilpy.ensure_installed`.

---

## Outils

| Script | Rôle | Réutilisable au-delà de l'éclipse |
|---|---|---|
| [`scripts/FusionHDR.py`](scripts/FusionHDR.py) | **Fusion HDR radiométrique** de poses bracketées (déjà alignées, linéaires) → FITS **32 bits** linéaire. Deux modes : *remplacement par seuil* (défaut) et *mélange pondéré SNR*. Échelle d'exposition lue dans l'en-tête ou estimée par recouvrement. | ✅ lunaire (clair de Terre, éclipse de Lune), planétaire (planète + satellites), **ciel profond** (cœur de M42…) |
| [`scripts/Corona.py`](scripts/Corona.py) | **Révélation de la couronne** : retrait du gradient radial (profil azimutal) + filtre **Larson–Sekanina** (radial/rotationnel) + accentuation (unsharp), avec sélection du centre et masque du disque lunaire. | ❌ spécifique à la couronne solaire |

> Les deux outils supposent des images **déjà alignées et calibrées** (l'alignement/calibration se fait en amont dans Siril). FusionHDR fusionne ≥ 2 poses ; Corona s'applique à **n'importe quelle** image de couronne (HDR ou non).

---

## Guides

Guide complet de bout en bout — **prise de vue → conversion → alignement/empilement par exposition → fusion HDR → couronne → finition** — incluant les **précautions de sécurité oculaire** (filtre retiré à C2 / remis à C3) :

- 🇫🇷 [Guide FR](docs/Guide_Eclipse_HDR_FR.md)
- 🇬🇧 [Guide US](docs/Guide_Eclipse_HDR_US.md)

---

## Installation

1. Copier `scripts/*.py` dans le dossier de scripts de Siril (ex. `~/siril/scripts/`).
2. Dans Siril : **Scripts → Rafraîchir**, puis lancer le script voulu.
3. Les dépendances Python manquantes sont installées automatiquement au premier lancement.

## Prérequis

- **Siril 1.3+**, **Python 3.10+** (via `sirilpy`)
- `FusionHDR` : `numpy`, `astropy`, `PyQt6`
- `Corona` : `numpy`, `astropy`, `scipy`, `PyQt6`

---

## Flux résumé

```
RAW bracketé → Siril (Conversion → FITS linéaires)
   → tri par exposition → align + empile par niveau (Siril/SolarAlign)
   → co-alignement des masters (Siril DFT)
   → FusionHDR (→ HDR 32 bits)        (facultatif si pas de bracketing)
   → Corona (retrait gradient radial + Larson–Sekanina)
   → Siril (Histogram / Asinh / GHS, couleur, export)
```

---

## Crédits

- Filtre **Larson–Sekanina** (gradient radial/rotationnel) ; méthode **Druckmüller** (NAFE/ACC) comme inspiration du rehaussement de couronne.
- Bibliothèques : numpy, astropy, scipy, PyQt6.

## Licence

[GPL-3.0-or-later](LICENSE). Chaque script porte l'en-tête `SPDX-License-Identifier: GPL-3.0-or-later`.
