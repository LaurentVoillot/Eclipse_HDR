# Eclipse_HDR

Outils Python pour [Siril](https://siril.org/) dédiés au **traitement HDR d'éclipse solaire** — et réutilisables pour le lunaire, le planétaire et le ciel profond.

*Python tools for Siril focused on solar-eclipse HDR processing — also usable for lunar, planetary and deep-sky imaging. The UI is in French; bilingual guides are included.*

Les scripts s'exécutent depuis **Siril → Scripts**, ouvrent une interface PyQt6 (thème sombre) et installent leurs dépendances au premier lancement via `sirilpy.ensure_installed`.

---

## Outils

| Script | Rôle | Réutilisable au-delà de l'éclipse |
|---|---|---|
| [`scripts/SirilJ_Align.py`](scripts/SirilJ_Align.py) | **Alignement, deux modes.** ☼ *Soleil* : recalage **multi-points** de la surface (corrélation de phase + champ de déformation local — corrige la turbulence, qu'une transformation globale ne peut pas suivre). ☾ *Éclipse* : recalage de poses de luminosités très différentes sur le **disque lunaire**, seul repère commun (détection adaptative : disque sombre enclos, ou limbe par RANSAC pour le diamant). Charge par défaut la séquence ouverte dans Siril. → 🇫🇷 [doc](docs/SirilJ_Align_FR.md) · 🇬🇧 [doc](docs/SirilJ_Align_US.md) | ✅ solaire haute résolution (mode ☼), lunaire |
| [`scripts/ConvertParVitesse.py`](scripts/ConvertParVitesse.py) | **Conversion RAW → FITS triés par vitesse.** Lit les EXIF sans rien convertir, affiche le plan de rangement, puis crée un sous-dossier par vitesse (1/160 s → `1-160`) avec les FITS dématriçés et la séquence. Les RAW ne sont ni copiés ni déplacés (liens temporaires). | ✅ toute série bracketée |
| [`scripts/FusionHDR_v1.2.py`](scripts/FusionHDR_v1.2.py) | **Variante expérimentale de FusionHDR.** Mode *Zones non saturées* : masque de saturation + marge + fondu **spatial** (le relais entre poses se fait loin de la saturation → pas de paliers). Masque lunaire continu, superposition des protubérances sans halo, coude doux en hautes lumières, calibration croisée. | ✅ idem FusionHDR |
| [`scripts/Corona_v1.1.py`](scripts/Corona_v1.1.py) | **Variante expérimentale de Corona.** Ajoute le *Détail tangentiel* (ACHF, flou le long des arcs — le seul exploitable sur champ bruité), le coude hautes lumières, la neutralisation de la couronne à Hα préservé, et le disque lunaire au niveau du fond. | ❌ couronne solaire |
| [`scripts/EclipseComposite.py`](scripts/EclipseComposite.py) | **Composition par couches** *(v0.1, esthétique assumée)* : Lune (clair de Terre), protubérances à sélectivité Hα, étoiles. | ❌ éclipse |
| [`scripts/FusionHDR.py`](scripts/FusionHDR.py) | **Fusion HDR radiométrique** de poses bracketées (déjà alignées, linéaires) → FITS **32 bits** linéaire. Deux modes : *remplacement par seuil* (défaut) et *mélange pondéré SNR*. Échelle d'exposition lue dans l'en-tête ou estimée par recouvrement. | ✅ lunaire (clair de Terre, éclipse de Lune), planétaire (planète + satellites), **ciel profond** (cœur de M42…) |
| [`scripts/Corona.py`](scripts/Corona.py) | **Révélation de la couronne** : retrait du gradient radial (profil azimutal) + filtre **Larson–Sekanina** (radial/rotationnel) + accentuation (unsharp), avec sélection du centre et masque du disque lunaire. | ❌ spécifique à la couronne solaire |

> **Chaîne recommandée :** ConvertParVitesse → SirilJ Align ☾ + `stack` (par vitesse) → SirilJ Align ☾ (co-alignement) → FusionHDR v1.2 → Corona v1.1 → GHS. Voir la [procédure complète](docs/Procedure_Complete_FR.md).
>
> Les scripts suffixés `_v1.1` / `_v1.2` sont des **variantes expérimentales** conservées à côté des originaux, qui restent inchangés — vous pouvez comparer les deux sur les mêmes données.

---

## Guides

**Procédure complète** — de la prise de vue à l'image finale : conversion, alignement et
empilement par exposition, fusion HDR, révélation de la couronne, finition. Inclut les
**précautions de sécurité oculaire** (filtre retiré à C2 / remis à C3), le dépannage, et
une annexe décrivant l'effet de chaque réglage de chaque écran.

- 🇫🇷 [Procédure complète FR](docs/Procedure_Complete_FR.md)
- 🇬🇧 [Full procedure EN](docs/Procedure_Complete_US.md)

Documentation par outil :

- **SirilJ Align** — 🇫🇷 [FR](docs/SirilJ_Align_FR.md) · 🇬🇧 [EN](docs/SirilJ_Align_US.md)

---

## Installation

1. Copier `scripts/*.py` dans le dossier de scripts de Siril (ex. `~/siril/scripts/`).
2. Dans Siril : **Scripts → Rafraîchir**, puis lancer le script voulu.
3. Les dépendances Python manquantes sont installées automatiquement au premier lancement.

## Prérequis

- **Siril 1.3+**, **Python 3.10+** (via `sirilpy`)
- `FusionHDR` : `numpy`, `astropy`, `PyQt6`
- `Corona` : `numpy`, `astropy`, `scipy`, `PyQt6`
- `SirilJ_Align` : `numpy`, `astropy`, `scipy`, `opencv-python`, `PyQt6`
- `ConvertParVitesse` : `exifread`, `PyQt6` (+ `exiftool` pour les RAW type CR3)

---

## Flux résumé

```
RAW bracketé → ConvertParVitesse (tri par vitesse + dématriçage)
   → SirilJ Align ☾ puis stack med/rej -nonorm -32b   (1 master par vitesse)
   → SirilJ Align ☾                                    (co-alignement des masters)
   → FusionHDR v1.2                                    (→ HDR 32 bits linéaire)
   → Corona v1.1                                       (couronne révélée)
   → EclipseComposite                                  (facultatif : Lune, étoiles)
   → Siril (GHS en mode luminance, saturation, export)
```

> ⚠️ Deux règles qui conditionnent tout le reste : **`set32bits` avant l'empilement**
> (sinon la couronne faible se terrasse), et **`-nonorm` à l'empilement** (toute
> normalisation casse la relation valeur ∝ pose × luminance dont dépend la fusion HDR).

---

## Crédits

- **RHEF** — Radial Histogram Equalizing Filter (Gilly & Cranmer, *Solar Physics* 2025)
- **FNRGF** — Fourier Normalizing Radial-Graded Filter (Druckmüllerová, Morgan & Habbal, *ApJ* 2011)
- **MGN** — Multi-Scale Gaussian Normalization (Morgan & Druckmüller, *Solar Physics* 2014)
- **ACHF** / détail tangentiel — d'après Druckmüller (*Contrib. Astron. Obs. Skalnaté Pleso* 2006)
- **Variance minimale** — Granados et al. 2010, reconstruction HDR optimale
- Filtre **Larson–Sekanina** (gradient radial/rotationnel), dans `Corona.py` d'origine
- Bibliothèques : numpy, astropy, scipy, opencv-python, exifread, PyQt6

## Licence

[GPL-3.0-or-later](LICENSE). Chaque script porte l'en-tête `SPDX-License-Identifier: GPL-3.0-or-later`.
