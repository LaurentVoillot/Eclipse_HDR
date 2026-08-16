# SirilJ Align — alignement Soleil & Éclipse

**Un seul outil, deux cibles.** Script Python pour Siril (`scripts/SirilJ_Align.py`, v1.0.0)
qui remplace et unifie deux outils antérieurs : l'alignement multi-points de la surface
solaire et le recalage des poses d'éclipse sur la Lune.

| Mode | Pour quoi | Repère utilisé |
|---|---|---|
| **☼ Soleil — surface** | séries à exposition constante (Hα, lumière blanche, granulation) | structure de surface, avec correction locale de la turbulence |
| **☾ Éclipse — Lune** | poses de luminosités très différentes (du diamant à la couronne) | disque lunaire, seul repère commun entre les niveaux |

Dans les deux modes, l'outil charge par défaut **la séquence ouverte dans Siril**.

---

## Pourquoi deux modes

Les deux problèmes sont opposés, d'où deux algorithmes.

**Soleil.** Toutes les images se ressemblent, mais la turbulence les déforme *localement* :
une transformation globale (translation, rotation) ne peut pas corriger une déformation qui
varie d'une zone à l'autre de l'image. Il faut mesurer un décalage en de nombreux points.

**Éclipse.** Les images ne se ressemblent pas du tout : entre une pose « diamant » et une
pose « couronne complète », le contenu diffère trop pour que la corrélation classique de
Siril fonctionne. Le **limbe lunaire** est le seul élément présent et identique partout.

---

## Prérequis

- Siril 1.3 ou plus récent, avec le support des scripts Python (`sirilpy`).
- Dépendances installées automatiquement : `numpy`, `astropy`, `scipy`, `opencv-python`, `PyQt6`.
- Les images doivent être des **FITS individuels**. Une séquence mono-fichier (FITSEQ/SER)
  est détectée et signalée : exportez-la en FITS individuels au préalable.

---

## ☼ Mode Soleil — surface (multi-points)

### Principe

1. **Translation globale** mesurée par corrélation de phase (précision sub-pixel).
2. **Grille de points d'ancrage** posée sur les zones réellement structurées de l'image de
   référence — granulation, protubérances, limbe. Le ciel vide est exclu automatiquement
   (un point n'est retenu que si son contenu dépasse le bruit).
3. **Décalage local** mesuré à chaque point, puis champ de déformation dense, lissé, avec
   rejet des points aberrants comblés par leurs voisins.
4. **Passes de raffinement** : chaque passe mesure le résidu et **compose les cartes** —
   l'image n'est donc rééchantillonnée **qu'une seule fois**, sans flou cumulé.

### Marche à suivre

1. Ouvrez la séquence dans Siril, puis lancez le script (il la détecte ; sinon **⟳ Rafraîchir**).
2. Choisissez l'**image de référence** : *Meilleure netteté* (défaut), *Première*, *Dernière*
   ou *Manuel*.
3. **Préparer la référence** → les points d'ancrage retenus s'affichent en bleu sur l'aperçu.
4. Ajustez si besoin (voir tableau), puis **▶ Démarrer l'alignement**.

### Réglages

| Réglage | Défaut | Rôle |
|---|---|---|
| Correction locale (multi-points) | coché | décoché = translation globale seule |
| Boîte (px) | 128 | taille du patch de mesure. Petit = suit un détail fin mais plus fragile ; grand = robuste mais lisse |
| Pas de grille (px) | 64 | espacement des points. Plus serré = champ plus fin, calcul plus long |
| Décalage local max (px) | 12 | au-delà, le point est rejeté (garde-fou anti-aberration) |
| Passes de raffinement | 2 | 2 suffit dans la plupart des cas ; l'image n'est rééchantillonnée qu'une fois |

Le nombre de points retenus s'affiche en direct — s'il est très faible, l'image est
probablement trop lisse : réduisez la boîte ou le pas.

### Sortie et intégration Siril

Les images alignées vont dans un **sous-dossier isolé** portant le préfixe choisi (défaut
`mp_`), à côté de la séquence. L'isolement est volontaire : la commande `convert` de Siril
mélangerait sinon les originaux et les images alignées dans la même séquence.

Options enchaînées automatiquement : **Convertir en séquence Siril**, **Empiler**
(`sum`, `med`, `rej 3 3`, `max`, `min`) et **Charger le résultat dans Siril**.

---

## ☾ Mode Éclipse — Lune (inter-poses)

### Principe

La détection du disque lunaire est **adaptative**, car une pose « couronne » et une pose
« diamant » n'offrent pas les mêmes indices :

- **Disque sombre enclos** par la couronne → centroïde de la composante sombre (précis,
  insensible à l'asymétrie de la couronne). Confiance 1,00.
- Sinon (diamant, poses très faibles) → **ajustement du limbe par RANSAC** puis affinage
  par moindres carrés (Kåsa), robuste même sur un arc partiel. Le pic le plus brillant est
  écrêté au préalable, sinon la lumière du diamant tire la détection vers elle.

Les poses sont ensuite recalées par translation pour faire coïncider les centres.

### Marche à suivre

1. Les poses de la séquence Siril courante se chargent automatiquement ; sinon
   **⟳ Séquence Siril courante**, ou **Ajouter** (FITS / Séquence / Dossier).
2. **↻ Détecter la Lune (toutes)**.
3. **Vérifiez les cercles**, en priorité sur les lignes de faible confiance
   (orange < 0,15, rouge < 0,07). C'est l'étape qui conditionne tout le reste.
4. **▶ Aligner sur la Lune**.

### Corriger un cercle à la main

| Geste | Effet |
|---|---|
| Clic-glissé sur l'aperçu | déplace le centre (validé au relâchement) |
| Maj + molette | ajuste le rayon |
| Molette | zoom · **Clic milieu glissé** : déplacer la vue |

Deux boutons font gagner beaucoup de temps :

- **Verrouiller au rayon médian** — impose le rayon médian des poses détectées à toutes.
- **⊕ Appliquer le centre affiché à toutes les poses** — la Lune bouge à peine à l'intérieur
  d'un bracket : réglez un centre parfait sur une pose nette, puis propagez-le.

### Sortie

Fichiers `aligned_00001.fit`, `aligned_00002.fit`… dans le dossier des poses (ou celui que
vous indiquez), plus la séquence **`aligned_.seq`** si l'option est cochée.

Deux points importants pour la suite du traitement :

- Les en-têtes sont **préservés**, `EXPTIME` compris — c'est ce dont FusionHDR a besoin.
- La sortie est divisée par un **facteur commun** pour tenir dans [0,1] : l'échelle
  **relative** entre poses est donc conservée, condition indispensable à une fusion HDR
  radiométriquement correcte.

Le fichier `aligned_.seq` est écrit directement plutôt que via `convert`, afin qu'il ne
référence **que** les images alignées, même si d'autres FITS traînent dans le dossier.

---

## Place dans la chaîne éclipse

```
ConvertParVitesse  →  SirilJ Align ☾  →  stack (par vitesse)  →  masters
                          ↓
                   SirilJ Align ☾  (co-alignement des masters)
                          ↓
                     FusionHDR  →  Corona  →  GHS
```

Le mode Éclipse sert **deux fois** :

1. **Avant chaque empilement**, sur les sous-poses d'un même niveau. C'est essentiel : si la
   Lune dérive pendant l'empilement, son bord est étalé et le master garde une morsure
   sombre définitive, qu'aucun traitement ultérieur ne peut rattraper.
2. **Sur les masters**, pour les co-aligner avant la fusion HDR.

> À l'empilement, utilisez `-nonorm` : toute normalisation casse la relation
> valeur ∝ pose × luminance dont dépend la fusion HDR.

Voir la [procédure complète](Procedure_Complete_FR.md).

---

## Dépannage

| Symptôme | Cause probable | Solution |
|---|---|---|
| « Aucune séquence chargée » | pas de séquence ouverte dans Siril | ouvrez-en une, puis **⟳ Rafraîchir** |
| « séquence mono-fichier » | FITSEQ ou SER | exportez en FITS individuels |
| Peu de points d'ancrage (Soleil) | image trop lisse ou trop bruitée | réduisez la boîte et le pas |
| « décalage global aberrant » | image très différente de la référence | changez de référence, ou écartez la pose |
| Cercle lunaire faux sur le diamant | détection tirée par le point brillant | corrigez au clic-glissé, puis **⊕ Appliquer à toutes** |
| Confiance faible sur plusieurs poses | couronne trop ténue pour cerner le disque | réglez un cercle sur une pose nette et propagez-le |

---

## Crédits

- Corrélation de phase multi-points — principe du *lucky imaging* solaire, à la manière d'AutoStakkert!
- Détection du limbe — ajustement de cercle RANSAC + moindres carrés de Kåsa
- Interface PyQt6, conventions VeraLux
- Licence GPL-3.0-or-later
