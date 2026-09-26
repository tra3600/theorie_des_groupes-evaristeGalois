"""Catalogue de groupes (et de faux groupes) prêts à l'emploi.

Toutes les permutations sont des tuples indexés à partir de 0 :
``p[i]`` est l'image de ``i``.  On les affiche en notation cyclique,
indexée à partir de 1, comme dans les livres : (1 2 3) envoie 1→2→3→1.
"""

from itertools import permutations
from math import gcd

from groupes import Groupe


# ----------------------------------------------------------------------
# Outils sur les permutations
# ----------------------------------------------------------------------
def composer(p, q):
    """p∘q : on applique d'abord q, puis p."""
    return tuple(p[i] for i in q)


def cycles(p):
    """Décomposition en cycles disjoints (indices à partir de 0, sans points fixes)."""
    vus, res = set(), []
    for depart in range(len(p)):
        if depart in vus or p[depart] == depart:
            continue
        cycle, i = [], depart
        while i not in vus:
            vus.add(i)
            cycle.append(i)
            i = p[i]
        res.append(cycle)
    return res


def notation_cyclique(p):
    """(0,2,1) -> '(2 3)' ; l'identité s'écrit 'id'."""
    cs = cycles(p)
    if not cs:
        return "id"
    return "".join("(" + " ".join(str(i + 1) for i in c) + ")" for c in cs)


def depuis_cycles(n, *cs):
    """Construit une permutation de {1..n} à partir de cycles écrits à partir de 1.

    >>> depuis_cycles(3, (1, 2, 3))
    (1, 2, 0)
    """
    p = list(range(n))
    for c in cs:
        for a, b in zip(c, c[1:] + c[:1]):
            p[a - 1] = b - 1
    return tuple(p)


def signature_permutation(p):
    """+1 pour une permutation paire, -1 pour une impaire."""
    return -1 if sum(len(c) - 1 for c in cycles(p)) % 2 else 1


# ----------------------------------------------------------------------
# Groupes classiques
# ----------------------------------------------------------------------
def cyclique(n):
    """ℤ/nℤ : les entiers modulo n avec l'addition."""
    return Groupe(range(n), lambda a, b: (a + b) % n, nom=f"Z/{n}Z")


def horloge():
    """ℤ/12ℤ déguisé en cadran d'horloge (0 s'affiche « 12h »)."""
    return Groupe(range(12), lambda a, b: (a + b) % 12, nom="Horloge",
                  etiquette=lambda a: f"{a or 12}h")


def inversibles_mod(n):
    """(ℤ/nℤ)* : les résidus premiers avec n, pour la multiplication."""
    elems = [a for a in range(1, n) if gcd(a, n) == 1] if n > 1 else [0]
    return Groupe(elems, lambda a, b: (a * b) % n, nom=f"(Z/{n}Z)*")


def symetrique(n):
    """S_n : toutes les façons de mélanger n objets."""
    return Groupe(permutations(range(n)), composer, nom=f"S{n}",
                  etiquette=notation_cyclique)


def alterne(n):
    """A_n : les permutations paires de S_n."""
    elems = [p for p in permutations(range(n)) if signature_permutation(p) == 1]
    return Groupe(elems, composer, nom=f"A{n}", etiquette=notation_cyclique)


def rotation(n, k):
    """Rotation du n-gone d'un angle 360·k/n (sommet i -> sommet i+k)."""
    return tuple((i + k) % n for i in range(n))


def reflexion(n, k):
    """Symétrie axiale du n-gone qui échange les sommets 0 et k."""
    return tuple((k - i) % n for i in range(n))


def nature_symetrie(p):
    """Décode une symétrie du polygone : ('rotation', k) ou ('miroir', k)."""
    n = len(p)
    k = p[0]
    return ("rotation", k) if p == rotation(n, k) else ("miroir", k)


def etiquette_polygone(p):
    genre, k = nature_symetrie(p)
    n = len(p)
    if genre == "rotation":
        return "id" if k == 0 else f"R{round(360 * k / n)}°"
    return f"M{k}"


def diedral(n):
    """D_n : les 2n symétries d'un polygone régulier à n côtés.

    Chaque symétrie est codée par la permutation des sommets qu'elle induit
    (sommet 0 en haut, numérotation dans le sens trigonométrique).
    """
    elems = [rotation(n, k) for k in range(n)] + [reflexion(n, k) for k in range(n)]
    return Groupe(elems, composer, nom=f"D{n}", etiquette=etiquette_polygone)


def klein():
    """Le groupe de Klein V4 = ℤ/2 × ℤ/2 (les symétries d'un rectangle)."""
    noms = {(0, 0): "e", (1, 0): "a", (0, 1): "b", (1, 1): "c"}
    return Groupe(noms, lambda x, y: ((x[0] + y[0]) % 2, (x[1] + y[1]) % 2),
                  nom="V4", etiquette=noms.get)


def _hamilton(p, q):
    a1, b1, c1, d1 = p
    a2, b2, c2, d2 = q
    return (a1 * a2 - b1 * b2 - c1 * c2 - d1 * d2,
            a1 * b2 + b1 * a2 + c1 * d2 - d1 * c2,
            a1 * c2 - b1 * d2 + c1 * a2 + d1 * b2,
            a1 * d2 + b1 * c2 - c1 * b2 + d1 * a2)


def quaternions():
    """Q8 = {±1, ±i, ±j, ±k} avec i² = j² = k² = ijk = -1 (Hamilton, 1843)."""
    unites = {"1": (1, 0, 0, 0), "i": (0, 1, 0, 0), "j": (0, 0, 1, 0), "k": (0, 0, 0, 1)}
    noms = {}
    for s, u in unites.items():
        noms[u] = s
        noms[tuple(-x for x in u)] = "-" + s
    return Groupe(noms, _hamilton, nom="Q8", etiquette=noms.get)


def melange_parfait(n_cartes, interieur=False):
    """Permutation d'un mélange « faro » parfait d'un paquet de n cartes (n pair).

    On coupe le paquet en deux moitiés égales puis on intercale les cartes
    une à une.  Mélange *extérieur* : la carte du dessus reste au-dessus.
    Mélange *intérieur* : la première carte de la 2e moitié passe dessus.
    Convention : ``nouveau_paquet[pos] = ancien_paquet[p[pos]]``.
    """
    if n_cartes % 2:
        raise ValueError("Il faut un nombre pair de cartes.")
    m = n_cartes // 2
    haut, bas = list(range(m)), list(range(m, n_cartes))
    premiers, seconds = (bas, haut) if interieur else (haut, bas)
    return tuple(c for paire in zip(premiers, seconds) for c in paire)


def appliquer_melange(paquet, p):
    return [paquet[i] for i in p]


def groupe_engendre_par_permutation(p, nom="<σ>"):
    """Le groupe cyclique {id, σ, σ², ...} sans énumérer tout S_n."""
    puissances, x = [tuple(range(len(p)))], p
    while x != puissances[0]:
        puissances.append(x)
        x = composer(p, x)
    return Groupe(puissances, composer, nom=nom, etiquette=notation_cyclique)


# ----------------------------------------------------------------------
# Faux groupes : ils échouent chacun à un axiome différent
# ----------------------------------------------------------------------
def faux_addition_bornee():
    """({-1, 0, 1}, +) : 1 + 1 = 2 sort de l'ensemble (pas de fermeture)."""
    return Groupe([-1, 0, 1], lambda a, b: a + b, nom="({-1,0,1}, +)")


def faux_soustraction(n=4):
    """(ℤ/nℤ, −) : interne, mais pas associative et sans vrai neutre."""
    return Groupe(range(n), lambda a, b: (a - b) % n, nom=f"(Z/{n}Z, -)")


def faux_multiplication_mod(n=6):
    """(ℤ/nℤ \\ {0}, ×) : n'est un groupe que si n est premier."""
    return Groupe(range(1, n), lambda a, b: (a * b) % n, nom=f"(Z/{n}Z\\{{0}}, ×)")


def faux_multiplication_avec_zero(n=5):
    """(ℤ/nℤ, ×) : 0 n'a jamais d'inverse."""
    return Groupe(range(n), lambda a, b: (a * b) % n, nom=f"(Z/{n}Z, ×)")


def faux_maximum():
    """({0,1,2,3}, max) : associative, neutre 0, mais max(3, x) ne vaut jamais 0."""
    return Groupe(range(4), max, nom="({0..3}, max)")
