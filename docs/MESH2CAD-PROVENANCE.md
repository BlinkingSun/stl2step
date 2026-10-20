# Provenance du fork Mesh2CAD

Initialisation : 20 octobre 2026.

| Élément | Référence |
| --- | --- |
| Dépôt amont | https://github.com/BlinkingSun/stl2step |
| Fork | https://github.com/Mycnc3d/stl2step |
| Commit amont de départ | `5ba491be6a77a4389b49255f12871dd4d4dfb13c` |
| Version moteur déclarée | 1.6.0 |
| Branche de préparation | `mesh2cad/bootstrap` |
| Prototype à migrer | https://github.com/Mycnc3d/MYCNC3D-Mesh2CAD |
| Commit du prototype inspecté | `d2de02ddf6927fadb938925636654233b4396ff1` |
| Version du prototype | 0.1.0.dev25 |

L'historique Git amont est conservé. Le ZIP initial était un support d'audit ;
le fork est fondé sur le dépôt Git, pas sur un nouvel historique reconstruit
depuis ce ZIP. Son SHA-256 était
`A40FA76E126BD4A64D6DD80C2065E73BDF22A84C7D453DA8CB3F860C62205B1C`.

## Licences et distribution

- Conserver le fichier LICENSE MIT et les mentions stl2step dans les sources.
- OCCT est une dépendance distincte, sous LGPL 2.1 avec exception additionnelle.
  Fournir les textes et sources correspondantes et vérifier les obligations
  applicables aux DLL effectivement distribuées.
- SolidOut est une application distincte dont les sources ne sont pas dans ce
  dépôt. Ne pas réutiliser son interface ou ses actifs propriétaires.
- Vérifier la licence du prototype et de ses dépendances avant de copier ses
  modules ; documenter leur provenance dans le commit de migration.
- Préparer un inventaire des dépendances du futur installateur.

## Synchronisation

`origin` désigne le fork Mycnc3d ; `upstream` désigne BlinkingSun/stl2step.
Les évolutions MYCNC3D sont préparées sur une branche dédiée. Les mises à jour
amont doivent être examinées, intégrées et testées avant livraison ; aucun
téléchargement automatique du moteur n'est nécessaire au fonctionnement.

## Validation de cette initialisation

Changements documentaires uniquement : aucun algorithme, test, fichier CMake
ou workflow de construction modifié. Aucune compilation ni conversion n'est
revendiquée pour cette initialisation. Les résultats annoncés par l'amont ne
constituent pas une qualification Windows MYCNC3D.
