"""Compter « à symétrie près » : le lemme de Burnside.

Deux coloriages sont « le même » si une symétrie du groupe G passe de l'un à l'autre.
Burnside : nombre de coloriages distincts = (1/|G|) · Σ_g (coloriages fixés par g)
                                            = (1/|G|) · Σ_g k^(nombre de cycles de g).
"""

from itertools import product

import catalogue as cat
from permgroup import points_du_cube, rotations_cube, action_sur, type_de_cycles


def nombre_coloriages(permutations, k):
    """Lemme de Burnside pour un groupe donné par la liste de TOUTES ses permutations."""
    total = sum(k ** len(type_de_cycles(p)) for p in permutations)
    assert total % len(permutations) == 0
    return total // len(permutations)


def indice_des_cycles(permutations):
    """Polynôme indicateur des cycles de Pólya : {type de cycles: nombre d'éléments}."""
    compte = {}
    for p in permutations:
        t = type_de_cycles(p)
        compte[t] = compte.get(t, 0) + 1
    return dict(sorted(compte.items(), reverse=True))


def orbites_force_brute(permutations, n_points, k):
    """Nombre d'orbites obtenu en listant les k^n coloriages (pour vérifier Burnside)."""
    vus, orbites = set(), 0
    for c in product(range(k), repeat=n_points):
        if c in vus:
            continue
        orbites += 1
        for p in permutations:
            vus.add(tuple(c[p[i]] for i in range(n_points)))
    return orbites


def colliers(n, k):
    """Colliers de n perles de k couleurs, à rotation près (groupe cyclique C_n)."""
    return nombre_coloriages([cat.rotation(n, j) for j in range(n)], k)


def bracelets(n, k):
    """Bracelets : à rotation ET retournement près (groupe diédral D_n)."""
    return nombre_coloriages(list(cat.diedral(n).elements), k)


def coloriages_du_cube(genre, k):
    """Coloriages de 'faces', 'sommets' ou 'aretes' d'un cube à rotation près (24 rotations)."""
    G, _ = rotations_cube()
    return nombre_coloriages([action_sur(genre, m) for m in G], k)


def groupe_des_rotations_du_cube():
    """Les 24 rotations du cube comme objet ``Groupe`` (éléments : matrices)."""
    from groupes import Groupe
    G, produit = rotations_cube()
    return Groupe(G, produit, nom="Rotations du cube", etiquette=lambda m: "".join(
        "+" if x > 0 else "-" if x < 0 else "0" for ligne in m for x in ligne))
