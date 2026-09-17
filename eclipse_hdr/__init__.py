# SPDX-License-Identifier: GPL-3.0-or-later
"""eclipse_hdr — traitement HDR d'éclipse solaire, du RAW à l'image finale.

Programme autonome : ne nécessite ni Siril ni aucune application externe.
Tous les RAW d'une session dans un dossier ; le pipeline en tire des masters
par vitesse, un HDR linéaire, la couronne révélée, puis les exports.

Chaque étape écrit un **TIFF 32 bits linéaire** avec ses métadonnées, peut
être lancée seule, interrompue, et reprise. Un export **16 bits** est
disponible à chaque étape pour Photoshop, DxO ou tout autre logiciel.

    from eclipse_hdr import Project, pipeline
    p = Project.open_or_create("~/totalite")
    pipeline.step_scan(p)
    pipeline.step_decode(p)
    pipeline.step_stack(p)

ou en ligne de commande :

    eclipse-hdr scan ~/totalite
    eclipse-hdr decode ~/totalite
    eclipse-hdr stack ~/totalite
"""

__version__ = "0.1.0"

from .project import Cancel, Project      # noqa: F401  (façade publique)

__all__ = ["Project", "Cancel", "__version__"]
