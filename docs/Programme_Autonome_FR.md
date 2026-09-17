# eclipse-hdr — programme autonome

Traitement HDR d'éclipse solaire **sans Siril** : tous les RAW dans un dossier,
un pipeline en étapes indépendantes, un TIFF 32 bits linéaire à chaque étape,
un export 16 bits à tout moment pour Photoshop, DxO ou Affinity.

> **État : phase 1.** Inventaire, décodage, empilement et export fonctionnent et
> sont validés sur RAW réels. Alignement, fusion HDR, couronne, composition et
> interface graphique restent à porter depuis les scripts Siril — voir
> [Ce qui reste à faire](#ce-qui-reste-à-faire).

---

## Installation

```bash
cd Eclipse_HDR
python3 -m venv .venv
.venv/bin/pip install -e .
```

Tout est installé : numpy, scipy, opencv, rawpy, tifffile, astropy, exifread.
Rien d'autre n'est requis — ni Siril, ni Photoshop, ni exiftool (celui-ci n'aide
que pour les EXIF des conteneurs récents type CR3 : `brew install exiftool`).

## Prise en main

```bash
.venv/bin/eclipse-hdr scan   ~/totalite
.venv/bin/eclipse-hdr decode ~/totalite
.venv/bin/eclipse-hdr stack  ~/totalite
.venv/bin/eclipse-hdr status ~/totalite
```

`~/totalite` est le **dossier projet**. Il contient les RAW, soit directement,
soit dans un sous-dossier `raw/`. Rien n'est écrit ailleurs, et les RAW ne sont
jamais modifiés, déplacés ni copiés.

## Le dossier projet

```
totalite/
  raw/                  ← vos RAW (ou directement dans totalite/)
  projet.json           ← manifeste : groupes, paramètres, état de chaque étape
  01_decode/1-160/…     ← linéaire + masque de saturation, rangé par vitesse
  03_stack/master_1-160.tif
  export/
```

Le manifeste est le **journal de bord** du traitement : il dit quelles vitesses
ont été trouvées, avec quels réglages chaque étape a tourné, quand, et combien
de temps. C'est lui qui rend la reprise possible.

**Reprise incrémentale.** Chaque étape retient une empreinte
`sha256(nom + paramètres + entrées)`. Relancer le pipeline ne recalcule que ce
dont l'empreinte a changé : retoucher un réglage de la couronne ne réempile
rien. Les fichiers déjà écrits sont sautés.

**Interruption.** Ctrl-C arrête après le fichier en cours (un second Ctrl-C
force). Toutes les écritures sont atomiques : une interruption ne laisse jamais
un fichier tronqué qui *aurait l'air* valide.

---

## Les étapes

| Commande | Rôle | Sortie |
|---|---|---|
| `scan` | lit les EXIF, groupe par vitesse | manifeste |
| `decode` | RAW → linéaire + masque de saturation | `01_decode/<vitesse>/` |
| `stack` | un master par vitesse | `03_stack/master_<vitesse>.tif` |
| `export` | TIFF 16 bits pour Photoshop/DxO | où vous voulez |
| `info` | métadonnées d'une image | terminal |
| `diff` | comparer deux images | terminal |
| `status` | état du projet | terminal |
| `run` | enchaîne les étapes disponibles | — |

### `scan` — inventaire

```bash
eclipse-hdr scan ~/totalite --tol 3
```

Lit le temps de pose de chaque RAW **sans rien convertir** : quelques secondes,
relançable librement pour voir le plan avant de s'engager. `--tol` est la
tolérance de regroupement en pour cent (3 % par défaut : fusionne les arrondis
APEX du type 1/1000 vs 1/1024, sans mélanger deux crans de ⅓ de diaphragme,
distants de 26 %).

Les vitesses deviennent des noms de dossier : 1/160 s → `1-160`, 2 s → `2s`,
1,3 s → `1s3`.

### `decode` — RAW vers linéaire

```bash
eclipse-hdr decode ~/totalite --speeds 1-160,2s --jobs 8
```

| Option | Défaut | Effet |
|---|---|---|
| `--speeds` | toutes | vitesses à traiter, séparées par des virgules |
| `--wb` | `daylight` | balance des blancs, la même pour toutes les images |
| `--jobs` | cœurs (max 8) | processus parallèles |
| `--demosaic` | `AHD` | `AHD`, `VNG`, `PPG`, `DHT`, `AAHD`, `LINEAR` |
| `--sat-frac` | `0.98` | seuil de saturation, en fraction de la plage utile |
| `--sat-dilate` | `3` | dilatation du masque, en px (portée du dématriçage) |
| `--half` | non | quart de résolution, pour un essai rapide |
| `--f32` | non | float32 au lieu de 16 bits (voir plus bas) |
| `--force` | non | refaire même si le fichier existe |

Compter environ **1,5 s par image et par cœur**. Quatre cents images sur huit
cœurs : à peu près une minute et demie.

#### `--wb` — fixer la balance des blancs

Une **seule** balance est appliquée à toutes les images de la session, et elle
est enregistrée dans le manifeste : les décodages suivants la réutilisent sans
qu'on ait à la redonner.

| Valeur | Balance retenue |
|---|---|
| *(rien)* ou `daylight` | « lumière du jour » de l'appareil — **défaut** |
| `camera` | telle que prise, lue sur le premier fichier |
| `neutral` | aucune : multiplicateurs égaux |
| `1.9722,0.9412,1.1376` | valeurs explicites, trois (RVB) ou quatre (RVBV) |
| `~/photos/_DSC0577.NEF` | balance de prise de vue de **ce** fichier |

Pour voir ce qu'un RAW contient avant de choisir :

```bash
eclipse-hdr info ~/totalite/raw/_DSC0577.NEF
```

```
pose : 1/640 s  (0.0015625 s)
balance lumière du jour  : 1.9722,0.9412,1.1376,0.9412   → --wb 1.9722,0.9412,1.1376
balance telle que prise  : 1.7930,1.0000,1.5352,1.0000   → --wb 1.7930,1.0000,1.5352
```

> **Seuls les rapports comptent.** libraw divise les quatre multiplicateurs par
> le plus petit avant de les appliquer : `2,2,2,2` donne exactement la même
> image que `1,1,1,1` (vérifié). Ce réglage change la **teinte**, jamais la
> luminosité.

Le défaut est la balance « lumière du jour » parce que c'est une **constante du
boîtier** — voir [L'invariant radiométrique](#linvariant-radiométrique) pour la
mesure qui l'impose.

**Changer d'avis en cours de route.** Si des images sont déjà décodées, un
`--wb` différent est **refusé** :

```
✗ la balance des blancs change (…) alors que 2 vitesse(s) sont déjà décodées.
  Mélanger deux balances dans un même projet casse l'échelle radiométrique
  dont dépend la fusion HDR.
  • pour tout redécoder avec la nouvelle balance : ajoutez --force
  • pour garder l'ancienne : relancez sans --wb
```

Avec `--force`, **toutes** les vitesses déjà décodées sont refaites, même celles
que `--speeds` ne demandait pas : un projet ne peut donc jamais se retrouver
avec deux balances mélangées. Une balance simplement proportionnelle à la
précédente n'est pas un conflit, puisqu'elle produit la même image.

### `stack` — un master par vitesse

```bash
eclipse-hdr stack ~/totalite
```

Sans `--method`, chaque groupe reçoit la méthode adaptée à son nombre de poses,
et le choix est journalisé :

| Poses | Méthode | Équivalent Siril |
|---|---|---|
| 2 – 3 | `median` | `stack med -nonorm` |
| 4 – 6 | `percentile` 0,2/0,1 | `stack rej p 0.2 0.1 -nonorm` |
| 7 et + | `winsorized` 3/3 | `stack rej w 3 3 -nonorm` |

| Option | Défaut | Effet |
|---|---|---|
| `--method` | auto | `median`, `mean`, `sum`, `min`, `max`, `percentile`, `sigma`, `winsorized` |
| `--low` / `--high` | selon méthode | seuils de rejet |
| `--missing` | `nan` | ce qui compte comme absence de donnée (voir plus bas) |
| `--band` | `256` | budget mémoire par bande, en Mo |

**Aucune normalisation n'est appliquée, et il n'y a aucun moyen d'en demander
une** : l'équivalent de `-nonorm` est le seul comportement possible. Le piège
qui casse la fusion HDR ne peut donc pas se refermer.

**Mémoire bornée.** Les images sont traitées par bandes horizontales : le pic
mémoire dépend de `--band`, pas du nombre de poses. Quarante images de 45 Mpx
tiennent dans 300 Mo au lieu de 21 Go.

### `export` — 16 bits pour Photoshop et DxO

```bash
eclipse-hdr export 03_stack/master_1-160.tif -o pour_photoshop.tif
eclipse-hdr export 05_hdr/hdr.tif --transform asinh --asinh-k 300
```

| `--transform` | Quand |
|---|---|
| `autostretch` *(défaut)* | usage courant — couleurs préservées (l'étirement porte sur la luminance) |
| `asinh` | étirement doux réglable par `--asinh-k`, garde le cœur non écrasé |
| `linear` | mise à l'échelle seule, si vous étirez ailleurs |

> ⚠ Le défaut n'est **pas** `linear`, à dessein. Un HDR linéaire couvre cinq
> décades ; ramené à 16 bits sans étirement, la couronne externe tombe sur
> quelques dizaines de niveaux et le fichier s'ouvre en rectangle noir.

---

## Les choix de conception, et pourquoi

Les trois affirmations ci-dessous ont été **mesurées** sur des RAW réels, pas
supposées. Les chiffres viennent d'un lot de vingt NEF.

### L'invariant radiométrique

Toute la fusion HDR repose sur *valeur ∝ pose × luminance*. Le décodage la
garantit par trois réglages : `gamma=(1,1)`, `no_auto_bright=True`, et des
multiplicateurs de balance des blancs **figés pour toute la session**.

La balance retenue est celle « lumière du jour » de l'appareil, une constante du
boîtier. Mesure : sur vingt fichiers, elle prend **une seule** valeur, tandis que
la balance *telle que prise* en prend **vingt différentes**. Utiliser cette
dernière donnerait à chaque image son propre facteur d'échelle.

Vérification directe : en divisant par deux le signal capteur, la sortie est
divisée par **0,4997** (écart 0,03 %), avec une dispersion de ±0,0004 sur
24,6 millions de pixels. Avec la luminosité automatique activée, le même essai
donne **×1,0006** — l'écart de pose est intégralement effacé.

### Le masque de saturation, lu côté capteur

Après balance des blancs et matrice couleur, un photosite saturé n'atterrit plus
à une valeur unique : chaque canal est multiplié différemment, puis mélangé.
Mesure : les pixels saturés au capteur ressortent étalés sur **0,254** en
sortie. Aucun seuil sur l'image finale ne peut donc les cerner proprement.

Le masque est lu **avant dématriçage**, par comparaison au niveau de saturation
du capteur, puis dilaté de trois pixels — un photosite saturé contamine ses
voisins interpolés (mesuré : 277 → 1855 pixels).

Le masque suit l'image jusqu'à la fusion, et l'empilement y écrit la **fraction
de poses saturées** en chaque pixel.

### Le stockage 16 bits du décodage

La sortie de libraw *est* du seize bits entier. La stocker en float32 la
rembourre de zéros : même information, deux fois le disque. Sur une session de
quatre cents images de 45 Mpx, cela fait **46 Go au lieu de 92**. L'aller-retour
est exact au bit près (vérifié).

Les étapes suivantes produisent du float32 : dès l'empilement, la moyenne de
plusieurs poses crée de vraies valeurs intermédiaires. `--f32` uniformise si
vous préférez.

---

## Validation contre Siril

Les mêmes quatre images (recadrées à 1200×800, soit 2 880 000 valeurs) empilées
par Siril 1.4.4 et par `eclipse-hdr`, à partir de données d'entrée **identiques
au bit près** :

| Méthode | Pixels divergents | Écart maximal |
|---|---|---|
| `median` | **0** | 0 — exact au bit près |
| `mean` | **0** | 3·10⁻⁸ (arrondi float32) |
| `percentile` 0,2/0,1 | **1 sur 2 880 000** | une valeur à 2·10⁻⁸ du seuil de rejet |
| `winsorized` 3/3 | divergent | voir ci-dessous |

Reproduire : `eclipse-hdr diff mon_master.tif master_siril.fit`.

### Deux divergences connues

**Les pixels nuls.** Siril exclut les pixels exactement nuls de ses empilements
par **moyenne** — mais pas de sa **médiane**. Sur `[0 ; 0,2179 ; 0,2313 ;
0,2459]` sa moyenne rend 0,2317 (la moyenne des trois non nuls) là où la
moyenne exacte vaut 0,1738. Cela concernait 4,3 % des pixels de l'essai.

La convention est défendable pour des poses recalées, où les bords laissés
vides par le décalage sont réellement sans donnée. Elle ne l'est pas pour des
poses brutes, où zéro est une mesure légitime (un photosite sous le niveau de
noir). Ici, une absence de donnée est donc marquée **NaN**, jamais zéro, et
`--missing zero` reproduit la convention de Siril — c'est avec elle que la
moyenne concorde exactement.

**L'écrêtage winsorisé.** Sur quatre poses, `stack rej w 3 3` de Siril ne rejette
*rien* : sa sortie est identique, bit pour bit, à sa moyenne simple. Celle d'ici
rejette, parce que le sigma winsorisé ne se laisse pas gonfler par la valeur
aberrante — c'est tout son intérêt.

La limite est chiffrable. Une valeur écartée de *d* parmi *n* poses porte
l'écart-type ordinaire à √(σ₀² + *d*²/*n*), et n'est rejetée à *k* σ que si
*d*²(1 − *k*²/*n*) > *k*²σ₀². Avec le *k* = 3 usuel le facteur s'annule dès que
**n ≤ 9** : sous dix poses, un écrêtage sigma ordinaire ne peut rejeter *aucune*
valeur isolée, si aberrante soit-elle. C'est la raison d'être de la
winsorisation, et l'explication de la règle « percentile sous sept poses ».

Cette divergence reste **à vérifier sur des poses d'éclipse** : quatre images de
scènes différentes ne sont pas un cas d'essai représentatif.

---

## Interopérabilité avec Siril

L'un ne remplace pas l'autre de force : les deux chaînes se croisent.

```bash
eclipse-hdr info  master.fit          # lit aussi les FITS
eclipse-hdr diff  a.tif b.fit         # compare TIFF et FITS
```

En Python, `imageio.from_fits` et `imageio.to_fits` importent et exportent des
FITS 32 bits avec `EXPTIME` et `MOONX`/`MOONY`/`MOONR` — de quoi injecter des
masters Siril à n'importe quelle étape, ou repartir dans Siril pour le GHS.

## Ce qui reste à faire

| Étape | État |
|---|---|
| Inventaire, décodage, empilement, export | **fait, validé sur RAW réels** |
| Alignement sur la Lune (☾) et multi-points (☼) | à porter depuis `SirilJ_Align.py` |
| Fusion HDR | à porter depuis `FusionHDR_v1.2.py` |
| Couronne (RHEF / FNRGF / MGN / ACHF) | à porter depuis `Corona_v1.1.py` |
| Composition par couches | à porter depuis `EclipseComposite.py` |
| Interface graphique PyQt6 | à écrire |

Les algorithmes existent déjà et sont indépendants de Siril : le portage
consiste à remplacer la couche d'entrées/sorties, pas à réécrire la science.

Le prétraitement (offsets, darks, flats) n'est pas prévu à ce stade, la chaîne
Siril actuelle ne l'utilisant pas non plus.

## Licence

[GPL-3.0-or-later](../LICENSE).
