# 🎓 Le laboratoire d'Évariste Galois

Un programme Python pour **découvrir la théorie des groupes en jouant** : des cas
illustrés (une horloge, un triangle en carton, un paquet de cartes, les quaternions…)
jusqu'au résultat qui a rendu Galois célèbre : **il n'existe pas de formule générale
pour résoudre les équations du 5e degré**.

Le programme ne se contente pas d'affirmer : chaque fait affiché est **calculé**
(tables, sous-groupes, isomorphismes, séries dérivées), et un quiz permet de s'entraîner.

![Les 6 symétries du triangle](figures/02_triangle_symetries.png)

## Installation

```sh
pip install -r requirements.txt     # matplotlib, networkx (+ pytest pour les tests)
```

## Utilisation

```sh
python groups-theory-EGalois.py                          # menu interactif
python groups-theory-EGalois.py --cas triangle           # un cas précis (répétable)
python groups-theory-EGalois.py --tout                   # tous les cas, figures à l'écran
python groups-theory-EGalois.py --tout --sauver figures  # figures enregistrées en PNG
python groups-theory-EGalois.py --tout --sans-graphique  # texte seulement
python groups-theory-EGalois.py --quiz 10 --graine 42    # le Défi de Galois (10 questions)
```

## Les 7 cas illustrés

| Cas | Idée | Ce que le programme montre |
|-----|------|----------------------------|
| `horloge` | 9h + 5h = 2h | Z/12Z, ses générateurs (1, 5, 7, 11), le cadran et son graphe de Cayley |
| `triangle` | reposer un triangle dans son trou | les 6 symétries de D3, leur table, la non-commutativité, D3 ≅ S3 |
| `carre` | les 8 symétries du carré | les 10 sous-groupes, le théorème de Lagrange, les classes, les sous-groupes distingués |
| `cartes` | le tour de magie du faro | 8 mélanges parfaits ramènent un jeu de 52 cartes dans l'ordre (52 en faro intérieur !) |
| `imposteurs` | des faux groupes | pour chacun, l'axiome violé **avec un contre-exemple** trouvé automatiquement |
| `quatre` | même taille, autre groupe | Z/4Z ≠ V4 ≅ (Z/8Z)*, D4 ≠ Q8 : l'« empreinte » des ordres des éléments |
| `galois` | l'équation du 5e degré | S3, S4 résolubles ; S5 bloqué sur A5, qui est simple (preuve par les classes de conjugaison) |

### L'horloge : Z/12Z

![Horloge](figures/01_horloge.png)

### Le carré : sous-groupes et graphe de Cayley

| Treillis des sous-groupes (en vert : distingués) | Graphe de Cayley |
|---|---|
| ![Sous-groupes de D4](figures/03_carre_sous_groupes.png) | ![Cayley D4](figures/03_carre_cayley.png) |

### Le mélange parfait de 52 cartes

Chaque ligne est le paquet après un mélange de plus : à la 8e, les couleurs sont revenues dans l'ordre.

![Faro](figures/04_cartes_52.png)

### Les imposteurs (extrait de la sortie)

```
▶ (Z/4Z, -)
  ✔ Fermeture      a·b reste toujours dans l'ensemble
  ✘ Associativité  (0·0)·1 = 3 mais 0·(0·1) = 1
  ✘ Neutre         aucun élément e tel que e·a = a·e = a
  ✘ Inverses       impossible sans élément neutre
  → Ce n'est PAS un groupe.
```

### Galois : pourquoi le degré 5 résiste

![Séries dérivées](figures/07_galois_series.png)

## Le Défi de Galois (quiz)

Des questions tirées au hasard dont la réponse est calculée par le programme :
heure sur l'horloge, inverse modulo p, ordre d'une permutation, composition de
symétries, Lagrange, mélanges de cartes… Indice avec `?` (½ point), séries 🔥, et
un titre final, d'« Apprenti·e calculateur·rice » à « Héritier·ère de Galois ».

```
Question 1/8 — Il est 9h. Quelle heure affichera l'horloge dans 17 heures ?
> 2
   ✔ Bravo !  9 + 17 = 26 ≡ 2 (mod 12) : on ne garde que le reste de la division par 12.
```

## Utiliser la bibliothèque

```python
from groupes import Groupe
import catalogue as cat

G = cat.diedral(4)                       # aussi : cyclique, horloge, symetrique, alterne,
print(G.verifier())                      #         inversibles_mod, klein, quaternions...
print(G.table_texte())
print(G.est_abelien(), G.contre_exemple_commutativite())
print(G.sous_groupes(), G.classes_de_conjugaison())
print(cat.symetrique(5).est_resoluble())  # False : merci Galois !

# Votre propre groupe (l'ancienne API is_group / identity / inverses marche toujours)
H = Groupe(range(5), lambda a, b: (a + b) % 5, nom="Z/5Z")
print(H.is_group(), H.identity, H.inverses)
```

## Organisation du code

| Fichier | Rôle |
|---------|------|
| `groups-theory-EGalois.py` | programme principal : menu, cas illustrés, options en ligne de commande |
| `groupes.py` | classe `Groupe` : axiomes avec contre-exemples, ordres, sous-groupes, classes, conjugaison, série dérivée, isomorphismes |
| `catalogue.py` | groupes prêts à l'emploi, permutations (notation cyclique), mélanges de cartes, faux groupes |
| `illustrations.py` | les figures matplotlib (tables de Cayley, graphes, polygones, horloge, cartes, treillis) |
| `quiz.py` | le Défi de Galois |
| `test_groupes.py` | tests (`python -m pytest -q`) |

Par rapport à la version d'origine, la table de Cayley est calculée une seule fois
(indices entiers) : vérifier les 1,7 million de triplets d'associativité de S5 prend
moins d'un dixième de seconde. Le graphe complet des 36 compositions, illisible, a été
remplacé par des graphes de Cayley construits sur des générateurs.
