# Construire Mesh2CAD sur le fork stl2step

## Décision

Le moteur stl2step est la fondation de la nouvelle application. Le prototype
Python/Tkinter fournit les fonctions d'interface à migrer, sans imposer une
réécriture du moteur C++ en Python. Windows x64 et le fonctionnement hors ligne
sont les premières cibles.

## Architecture cible

```text
Interface MYCNC3D : visualisation, unités, mesures, sections, SpaceMouse
    -> tâche de conversion : entrée figée, annulation, fichiers temporaires
    -> moteur stl2step local : C++ / OpenCASCADE
    -> validation indépendante du STEP et rapport
    -> aperçu du résultat puis remplacement atomique de la destination
```

Le moteur reste à la racine avec ses interfaces actuelles. La future interface
sera placée sous `desktop/`, séparée des algorithmes natifs. Un processus moteur
par travail évite de partager la configuration globale OCCT avec l'interface.
Les unités sont explicites ; le STEP est en millimètres.

Deux modes seront présentés : conserver les facettes (Verbatim) ou reconstruire
les surfaces mécaniques (TrueForm). Le mode analytique ne restitue pas un arbre
de conception paramétrique et ne garantit pas une reconstruction de scans.

## Étapes et critères d'acceptation

### 1. Qualifier le moteur Windows

- Reprendre la chaîne CMake/MSVC existante et figer OCCT. Le workflow Windows
  amont référence OCCT 8.0.1 et vcpkg
  `f89a4a1da4e3176a8d1a14c1825b9b2f98e48843`.
- Construire et exécuter les tests avec cette configuration ; enregistrer les
  versions et les résultats, sans modifier les seuils pour masquer un échec.
- Comparer cube, cylindre, platine percée, pièce à niveaux, multicorps, maillage
  ouvert et cas ambigu hexagone/cercle. Ajouter des pièces réelles autorisées.
- Contrôler les surfaces du STEP relu, pas uniquement les surfaces reconnues.

Livrable : binaire Windows qualifié et rapport reproductible ; les performances
et la mémoire sont mesurées sur les mêmes fichiers et la même machine.

### 2. Migrer une première interface utilisable

- Récupérer du prototype le viewport, la caméra, les unités et la SpaceMouse
  après vérification des licences ; conserver les tests utiles.
- Brancher le lancement du moteur sur la gestion de travaux existante.
- Utiliser une liste d'arguments, jamais une commande shell construite à partir
  des chemins. Consommer uniquement la ligne `RESULT ` comme contrat JSON.
- Préserver les journaux, version du moteur, options et empreinte de l'entrée.
- Prévoir annulation, erreurs, chemins Unicode, espaces et destination existante.
- Pour OBJ/PLY, définir une conversion intermédiaire avec contrôle des unités
  et de la précision ; le moteur accepte directement le STL uniquement.

Livrable : ouvrir un STL, choisir le mode, convertir, annuler et afficher le
rapport sans bloquer l'interface ni perdre une destination existante.

### 3. Valider avant de livrer le STEP

`ok=true` signifie qu'un fichier a été produit, pas qu'il est prêt à usiner.
Les sorties 0 et 2, les avertissements et enveloppes ouvertes restent distincts.
Ne pas publier le fichier temporaire sur le seul critère du code de sortie.

La validation doit vérifier : relecture, B-Rep valide, fermeture, corps attendus,
dimensions, volumes et surfaces réellement présentes. Les écarts locaux doivent
être mesurés séparément ; un volume identique ne prouve pas la fidélité locale.
`volumeDeltaPct` mesure le transfert B-Rep vers STEP, pas la fidélité au maillage.
Un contrôle non effectué doit être indiqué comme tel, jamais présenté comme
réussi. Le moteur amont saute notamment un contrôle global au-delà de 50 000
faces. Conserver sa vérification tant que le validateur indépendant n'est pas
opérationnel.

Livrable : statuts explicites valide, à examiner et échec ; remplacement atomique
par le processus parent seulement après validation et absence d'annulation.

### 4. Prioriser les pièces prismatiques et le DXF

Qualifier les profils lignes/arcs, trous et extrusions par niveau. Réutiliser
`refit_prism*`, `refit_profile*` et `dxf_export` via le moteur avant d'en extraire
des composants. Migrer les sections et mesures du prototype selon les besoins.

Livrable : cas de référence avec rayons, diamètres et dimensions contrôlés.

### 5. Distribuer hors ligne

Créer un paquet autonome contenant interface, moteur, DLL et notices. Tester
sur Windows 10/11 x64 sans Python préinstallé ni accès Internet. La construction
peut télécharger des dépendances ; l'application livrée doit fonctionner sans
connexion, compte, activation distante ou mise à jour obligatoire.

## Politique de développement

Garder de petits commits séparant interface, corrections moteur et mise à jour
amont. Respecter les contrôles de contribution d'AGENTS.md pour les changements
affectant la compilation. Le présent démarrage ne modifie que la documentation.
Le prototype d'origine reste conservé ; sa migration n'implique ni suppression
ni écrasement de ses modifications locales.
