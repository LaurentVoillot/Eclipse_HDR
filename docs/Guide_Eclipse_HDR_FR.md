# Guide — Traitement HDR d'une éclipse solaire

**Du RAW à la couronne finie dans Siril.**
Chaîne complète : prise de vue → conversion → alignement/empilement par exposition → fusion HDR → rehaussement de couronne → finition.

---

## ⚠ Sécurité — vue & matériel (à lire en premier)

**En dehors de la totalité, le Soleil reste dangereux pour les yeux et le capteur.** La photosphère est visible pendant toutes les phases partielles (avant C2, après C3) : la regarder — ou la viser — sans protection cause des **brûlures rétiniennes irréversibles et indolores**, et peut détruire le capteur ou provoquer un incendie.

| Moment | Yeux | Optique / appareil |
|---|---|---|
| **Phases partielles** (avant C2, après C3) | **Lunettes éclipse certifiées ISO 12312-2** uniquement. Jamais à l'œil nu, ni dans un viseur / des jumelles / un télescope sans filtre. | **Filtre solaire certifié** (ex. Baader AstroSolar) **devant** l'objectif, en permanence. |
| **Totalité** (entre C2 et C3, **seulement**) | Observation à l'œil nu **autorisée** — c'est le seul moment. | Photo **sans filtre** autorisée. |

- **À C2 / fin du diamant** : retirez le filtre **seulement** une fois les grains de Baily et le diamant disparus, totalité pleinement installée.
- **À C3 / réapparition du diamant** : la photosphère revient **brutalement**. **Remettez le filtre et cessez toute observation directe à l'instant où le diamant réapparaît** — n'attendez pas : une fraction de seconde suffit à blesser.
- Connaissez vos **horaires de contact (C2, C3)** précis, gardez le filtre **prêt à remettre**, et prévenez toute personne présente. En cas de doute : **filtre en place**.

> Les lunettes éclipse protègent **les yeux**, pas l'appareil ; le filtre d'objectif protège **l'appareil**, pas vos yeux si vous regardez ailleurs. Les deux sont nécessaires pendant les phases partielles.

---

## Outils de la chaîne

| Étape | Outil |
|---|---|
| Conversion RAW → FITS | **Siril** (onglet *Conversion*) |
| Alignement + empilement par niveau | **Siril** (*Registration* + *Stacking*) ou **SolarAlign** |
| Co-alignement des masters | **Siril** (registration DFT) |
| Fusion HDR 32 bits | **FusionHDR** |
| Révélation de la couronne | **Corona** |
| Courbes / finition | **Siril** (Histogram / Asinh / GHS) |

> Pré-requis logiciels : Siril 1.3+, et les scripts `FusionHDR.py` / `Corona.py` dans le dossier de scripts de Siril.

---

## Phase 0 — Prise de vue (le jour J)

La qualité du résultat se joue ici.

- **Monture équatoriale** en suivi, bien mise en station → le Soleil reste cadré, l'alignement est presque pure translation (pas de rotation de champ).
- **Mise au point manuelle** verrouillée (faites-la avant sur le limbe lunaire ou une étoile brillante), puis **scotchez la bague**. Ne touchez plus au zoom/focus.
- **Exposition manuelle** : ISO fixe (100–400), diaphragme fixe, on ne fait varier que la **vitesse**.
- **Bracketing large** : de ~**1/4000 s à ~2–4 s**, par pas de ~1 diaphragme → **~12–15 niveaux**. La couronne couvre un énorme rapport de luminosité ; il faut toute la gamme.
- **RAW** obligatoire (données linéaires). Rafales continues, intervallomètre ou bracketing auto du boîtier, **cadrage identique**.
- **Filtre solaire** : retirez-le **à C2** (totalité pleinement installée, diamant disparu) et **remettez-le à C3 dès la réapparition du diamant** — répétez ce geste à blanc avant le jour J. Voir la section **⚠ Sécurité** : c'est vital pour vos yeux **et** le capteur.
- Pensez aux instants spécifiques : **grains de Baily / diamant** (poses très courtes) et **protubérances** (poses courtes), juste après C2 et avant C3.
- Anti-vibration : retardateur, obturateur électronique ou relevage du miroir.

> ⚠ **La Lune se déplace devant la couronne** (~0,5″/s, soit ~2′ sur 4 min). Voir « Pièges » : cela conditionne le choix d'alignement.

---

## Phase 1 — Conversion des RAW en FITS (Siril)

1. Onglet **Conversion**.
2. Glissez tous les RAW (ou par groupe — voir Phase 2).
3. Sortie **FITS**, cochez **Debayer** si appareil **couleur** (OSC/DSLR).
4. **Convertir**. Siril produit des FITS linéaires (les `EXPTIME`/EXIF sont conservés — utile pour FusionHDR).

> Gardez les données **linéaires**, n'étirez rien à ce stade.

---

## Phase 2 — Tri par niveau d'exposition

La fusion HDR combine **un master par luminosité**. Il faut donc regrouper les poses par temps de pose.

- Classez les RAW dans des **dossiers par exposition** : `1_2000/`, `1_500/`, `1_125/`, … `2s/` (d'après la vitesse EXIF).
- Convertissez **chaque dossier** en sa propre séquence FITS (Phase 1 par dossier).

Résultat : une séquence FITS par niveau d'exposition.

---

## Phase 3 — Alignement + empilement de chaque niveau

Pour **chaque** séquence (chaque luminosité) :

**Option A — Siril natif (recommandé)**
1. Onglet **Registration** → méthode **« Image Pattern Alignment (DFT) »** (corrélation de phase, translation). Elle accroche le motif dominant à fort contraste : **disque lunaire + couronne interne**. (La méthode « Global Star Alignment » ne marche pas : pas d'étoiles.)
2. **Register**.
3. Onglet **Stacking** → **Moyenne avec rejet** (Winsorized/Sigma, ex. 3/3) pour éliminer avions, satellites, rafales nuageuses.
4. **Stack** → un **master** pour ce niveau (ex. `master_1_500.fit`).

**Option B — script SolarAlign**
Chargez la séquence du niveau dans Siril, lancez **SolarAlign** (StackReg *Corps rigide* ou ECC *Euclidien*), cochez *Convertir* + *Empiler* → il aligne et empile en une fois. Pratique si la registration native peine.

> **Astuce poses longues vs courtes** :
> - **Poses longues** (couronne externe, faible) → empilez-en beaucoup (réduit le bruit) ; le bord lunaire surexposé sera masqué de toute façon.
> - **Poses courtes** (protubérances, couronne interne) → **fenêtre temporelle courte** ou meilleures images uniques (les protubérances évoluent, la Lune dérive).

---

## Phase 4 — Co-alignement des masters

Les masters doivent **se superposer** avant la fusion HDR (sinon artefacts de bord).

1. Placez tous les masters dans un même dossier, convertissez-les en **une séquence**.
2. **Registration** → **DFT** (le disque lunaire, commun à toutes les poses, sert d'ancrage malgré les luminosités différentes).
3. Exportez/enregistrez les images **recalées** (préfixe `r_`).

> Si la monture a parfaitement suivi avec un cadrage strictement identique, les masters peuvent déjà être alignés ; vérifiez, sinon faites cette phase.

---

## Phase 5 — Fusion HDR (FusionHDR)

> **Pas de bracketing ?** Cette phase est **facultative**. Si vous n'avez qu'une seule exposition (un seul master empilé), **sautez la Phase 5** et passez directement à **Corona** (Phase 6).

1. **Scripts → FusionHDR**.
2. **Ajouter…** les masters recalés (un par exposition).
3. Vérifiez la colonne **Expo (s)** (lue dans l'en-tête). Si absente/identique → **Estimer les expositions (recouvrement)**.
4. **Mode** : *Remplacement par seuil* (défaut, le plus prévisible).
5. **Pleine échelle** : *Auto*. **Seuil saturation** ~0,95.
6. **Soustraire le fond** : **décoché** pour l'éclipse (ciel ~noir).
7. **Normaliser la sortie** : coché.
8. **Fusionner** → aperçu. **Enregistrer + charger dans Siril** → un **FITS 32 bits linéaire** HDR.

---

## Phase 6 — Révélation de la couronne (Corona)

Sur le HDR chargé dans Siril :

1. **Scripts → Corona**.
2. **Centre** : cliquez **Auto**, puis affinez au **clic** sur le centre du disque. Réglez le **rayon lunaire** pour que le **cercle bleu** colle au bord de la Lune.
3. **Retrait du gradient radial** : **activé** (c'est l'étape clé — elle aplatit la décroissance et fait ressortir les streamers). Lissage ~8.
4. **Appliquer** → la couronne apparaît déjà structurée.
5. Optionnel **Larson–Sekanina** : Δ angle ~3–8°, Δ radial 0, intensité ~1 → crispe les **streamers** (structures radiales). Augmentez prudemment (amplifie le bruit).
6. Optionnel **Unsharp** : sigma 2–4, intensité ~1 → détail fin.
7. **Noircir le disque lunaire** : coché.
8. **Enregistrer + charger dans Siril**.

> Itérez : ajustez centre / rayon / paramètres et **Appliquer** à nouveau jusqu'au rendu voulu.

---

## Phase 7 — Finition dans Siril (les courbes)

Le résultat de Corona est aplati : place à l'esthétique.

- **Histogram Transformation** (étirement) ou **Asinh** pour révéler la dynamique douce de la couronne.
- **Generalized Hyperbolic Stretch (GHS)** pour un contrôle fin du contraste local sans cramer le centre.
- **Couleur** : balance des blancs / *Photometric Color Calibration* si pertinent ; saturation modérée (la couronne est légèrement structurée en couleur).
- **Réduction de bruit** légère sur les zones faibles si besoin.
- Export final : **TIFF 16 bits** ou **PNG**.

---

## Pièges & variantes

- **Lune vs couronne** : la couronne est fixe (Soleil), la Lune dérive. Sur une totalité longue, traitez **deux versions** — une **alignée couronne** (poses longues, couronne externe) et une **alignée Lune** (poses courtes, protubérances/couronne interne) — puis **compositez** (la couronne + la Lune/protubérances rephotographiées au bon instant). Pour une première passe ou une totalité courte, l'alignement DFT sur le disque lunaire suffit.
- **Protubérances & chromosphère** : viennent des poses **courtes** d'un instant précis (juste après C2 / avant C3), **non empilables** sur toute la totalité. Ajoutez-les en calque/composite.
- **Disque lunaire** : souvent remplacé par un disque noir propre, ou par une **lune en clair de Terre** photographiée séparément (pose longue).
- **Couleur** : FusionHDR et Corona gèrent le couleur (modèle radial sur la luminance, application par canal) ; gardez la balance en fin de chaîne.
- **Sans HDR (une seule exposition)** : la chaîne marche aussi **sans bracketing**. Empilez vos poses (Phase 3), puis appliquez **Corona** directement sur le master en **sautant la Phase 5** (FusionHDR). Corona révèle la couronne sur **n'importe quelle** image — HDR ou non — du moment qu'elle n'est pas trop saturée au centre. (FusionHDR, lui, requiert au moins deux poses de luminosités différentes.)
- **Mémoire** : FusionHDR fusionne **une image à la fois** (tient sur de gros fichiers) ; Corona charge l'image entière.

---

## Récapitulatif express

```
RAW (bracketing, RAW, monture en suivi, filtre retiré à C2)
   │  Siril : Conversion (+ Debayer si couleur) → FITS linéaires
   ▼
Tri par exposition (dossiers) 
   │  par niveau : Registration DFT + Stacking (rejet) → 1 master/expo
   ▼  (ou SolarAlign : align + empile)
Co-alignement des masters (Siril DFT) → masters recalés
   │
   ▼  FusionHDR : remplacement par seuil, échelle auto → HDR 32 bits   (facultatif — sauter si pas de bracketing)
   │
   ▼  Corona : centre + retrait gradient radial (+ LS + unsharp) → couronne révélée
   │
   ▼  Siril : Histogram / Asinh / GHS, couleur, export TIFF/PNG
RÉSULTAT
```

---

*Cette même chaîne FusionHDR + Corona se réutilise au-delà de l'éclipse : FusionHDR pour le lunaire (clair de Terre, éclipse de Lune), le planétaire (planète + satellites) et le ciel profond (cœur de M42…) ; Corona reste spécifique à la couronne solaire.*
