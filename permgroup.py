"""Grands groupes de permutations : ordre, appartenance, cube de Rubik.

La classe ``Groupe`` stocke toute la table de Cayley : impossible pour le Rubik's Cube
(43 252 003 274 489 856 000 éléments). ``GroupePermutations`` ne manipule que quelques
générateurs et la *chaîne de stabilisateurs* de Schreier–Sims : on lit alors l'ordre du groupe
comme un produit de tailles d'orbites, sans jamais énumérer ses éléments.

L'algorithme est la version aléatoire de Schreier–Sims : exact quand il répond « trouvé »,
et la probabilité d'un ordre trop petit est inférieure à 2⁻³⁰ avec les réglages par défaut.
``enumerer`` donne un calcul exact (mais exponentiel) pour recouper sur de petits groupes.
"""

import random
from math import gcd

Perm = tuple


def identite(n):
    return tuple(range(n))


def compose(p, q):
    """p∘q : on applique d'abord q, puis p (même convention que ``catalogue.composer``)."""
    return tuple(p[i] for i in q)


def inverse(p):
    inv = [0] * len(p)
    for i, x in enumerate(p):
        inv[x] = i
    return tuple(inv)


def ordre_permutation(p):
    """Plus petit k tel que p^k = id : le ppcm des longueurs des cycles."""
    vus, o = set(), 1
    for d in range(len(p)):
        if d in vus:
            continue
        longueur, i = 0, d
        while i not in vus:
            vus.add(i)
            i = p[i]
            longueur += 1
        o = o * longueur // gcd(o, longueur)
    return o


def type_de_cycles(p):
    """Longueurs des cycles, y compris les points fixes, triées de la plus grande à la plus petite."""
    vus, res = set(), []
    for d in range(len(p)):
        if d in vus:
            continue
        longueur, i = 0, d
        while i not in vus:
            vus.add(i)
            i = p[i]
            longueur += 1
        res.append(longueur)
    return tuple(sorted(res, reverse=True))


class GroupePermutations:
    """Groupe engendré par des permutations d'un même ensemble {0, ..., n-1}."""

    def __init__(self, generateurs, graine=0):
        gens = [tuple(g) for g in generateurs]
        if not gens:
            raise ValueError("Il faut au moins un générateur.")
        self.degre = len(gens[0])
        if any(sorted(g) != list(range(self.degre)) for g in gens):
            raise ValueError("Les générateurs doivent être des permutations de {0..n-1}.")
        self.generateurs = gens
        self._id = identite(self.degre)
        self._rng = random.Random(graine)
        self._base = []
        self._fortes = []        # générateurs forts
        self._orbites = []       # par niveau : {point: transversal u tel que u(base[i]) = point}
        self._construit = False

    # ------------------------------------------------------------------
    # Chaîne de stabilisateurs
    # ------------------------------------------------------------------
    def _recalculer(self):
        self._orbites = []
        for i, b in enumerate(self._base):
            fixant = [s for s in self._fortes if all(s[x] == x for x in self._base[:i])]
            orbite = {b: self._id}
            a_voir = [b]
            while a_voir:
                x = a_voir.pop()
                for s in fixant:
                    y = s[x]
                    if y not in orbite:
                        orbite[y] = compose(s, orbite[x])
                        a_voir.append(y)
            self._orbites.append(orbite)

    def _filtrer(self, g):
        """Fait passer g dans la chaîne : (résidu, niveau). Résidu = id et niveau = len(base) ssi g ∈ G."""
        for i, b in enumerate(self._base):
            beta = g[b]
            if beta not in self._orbites[i]:
                return g, i
            g = compose(inverse(self._orbites[i][beta]), g)
        return g, len(self._base)

    def _ajouter(self, residu):
        if all(residu[x] == x for x in self._base):
            self._base.append(next(x for x in range(self.degre) if residu[x] != x))
        self._fortes.append(residu)
        self._recalculer()

    def _construire(self, succes_requis=30):
        if self._construit:
            return
        self._base, self._fortes = [], []
        for g in self.generateurs:
            if g != self._id:
                if not self._base or self._filtrer(g) != (self._id, len(self._base)):
                    self._ajouter(self._filtrer(g)[0] if self._base else g)
        if not self._base:
            self._construit = True
            return
        pool = list(self.generateurs) * 3 + [self._id]
        succes, essais = 0, 0
        while succes < succes_requis:
            essais += 1
            i, j = self._rng.sample(range(len(pool)), 2)
            pool[i] = compose(pool[i], pool[j])           # remplacement de produit
            r, niveau = self._filtrer(pool[i])
            if r != self._id:
                self._ajouter(r)
                succes = 0
            else:
                succes += 1
        self._construit = True

    # ------------------------------------------------------------------
    # Interface
    # ------------------------------------------------------------------
    def ordre(self):
        """Nombre d'éléments du groupe (produit des tailles d'orbites de la chaîne)."""
        self._construire()
        n = 1
        for o in self._orbites:
            n *= len(o)
        return n

    @property
    def base(self):
        self._construire()
        return list(self._base)

    def tailles_orbites(self):
        """Tailles des orbites de la chaîne : leur produit est l'ordre du groupe."""
        self._construire()
        return [len(o) for o in self._orbites]

    def contient(self, p):
        """p appartient-il au groupe ?"""
        self._construire()
        p = tuple(p)
        if len(p) != self.degre:
            return False
        r, niveau = self._filtrer(p)
        return r == self._id and niveau == len(self._base)

    def orbite(self, point):
        """Orbite d'un point sous l'action du groupe."""
        vus, a_voir = {point}, [point]
        while a_voir:
            x = a_voir.pop()
            for g in self.generateurs:
                if g[x] not in vus:
                    vus.add(g[x])
                    a_voir.append(g[x])
        return sorted(vus)

    def est_transitif(self):
        return len(self.orbite(0)) == self.degre

    def element_aleatoire(self):
        """Un élément (quasi) uniforme, obtenu en recomposant des éléments de la chaîne."""
        self._construire()
        g = self._id
        for orbite in reversed(self._orbites):
            g = compose(g, self._rng.choice(list(orbite.values())))
        return g

    def enumerer(self, limite=300_000):
        """Tous les éléments, par parcours en largeur (exact, mais seulement pour les petits groupes)."""
        vus = {self._id}
        a_voir = [self._id]
        while a_voir:
            x = a_voir.pop()
            for g in self.generateurs:
                y = compose(g, x)
                if y not in vus:
                    vus.add(y)
                    a_voir.append(y)
                    if len(vus) > limite:
                        raise ValueError(f"plus de {limite} éléments : utiliser ordre()")
        return vus

    def distribution_des_ordres(self, echantillon=500):
        """Ordres de ``echantillon`` éléments tirés au hasard : {ordre: effectif}."""
        compte = {}
        for _ in range(echantillon):
            o = ordre_permutation(self.element_aleatoire())
            compte[o] = compte.get(o, 0) + 1
        return dict(sorted(compte.items()))


# ----------------------------------------------------------------------
# Cubes de Rubik, construits par la géométrie
# ----------------------------------------------------------------------
def _rotation(axe, signe):
    """Quart de tour (de matrice entière) autour de l'axe 0/1/2, dans un sens ou dans l'autre."""
    if signe not in (1, -1):
        raise ValueError("signe = ±1")
    if axe == 0:
        return lambda v: (v[0], -signe * v[2], signe * v[1])
    if axe == 1:
        return lambda v: (signe * v[2], v[1], -signe * v[0])
    return lambda v: (-signe * v[1], signe * v[0], v[2])


def autocollants(n=3):
    """Les autocollants d'un cube n×n×n (n = 2 ou 3) : couples (position du cubie, normale)."""
    pos = [-1, 1] if n == 2 else [-1, 0, 1]
    res = []
    for x in pos:
        for y in pos:
            for z in pos:
                c = (x, y, z)
                if c == (0, 0, 0):
                    continue
                for axe in range(3):
                    if abs(c[axe]) == 1:
                        normale = tuple(c[axe] if k == axe else 0 for k in range(3))
                        res.append((c, normale))
    return res


def mouvement_cube(face, n=3):
    """Permutation des autocollants pour un quart de tour de la face R, L, U, D, F ou B
    (dans le sens des aiguilles d'une montre vue de l'extérieur)."""
    axe = "RLUDFB".index(face) // 2
    cote = 1 if face in "RUF" else -1
    sens = -cote            # horaire vu de l'extérieur = rotation de signe -côté autour de l'axe
    tourner = _rotation(axe, sens)
    st = autocollants(n)
    index = {s: i for i, s in enumerate(st)}
    p = list(range(len(st)))
    for i, (c, nor) in enumerate(st):
        if c[axe] == cote:
            p[i] = index[(tourner(c), tourner(nor))]
    return tuple(p)


def cube_rubik(n=3, faces=None):
    """Groupe du cube n×n×n engendré par des quarts de tour.

    n = 3 : les six faces. n = 2 : R, U, F seulement (on garde un coin fixe, sinon on
    compterait aussi les 24 orientations du cube entier)."""
    faces = faces or ("RLUDFB" if n == 3 else "RUF")
    return GroupePermutations([mouvement_cube(f, n) for f in faces])


def rotations_cube():
    """Les 24 rotations du cube, comme matrices 3×3 entières (tuples de lignes)."""
    def mat(f):
        return tuple(tuple(f(tuple(1 if i == j else 0 for i in range(3)))[k] for j in range(3))
                     for k in range(3))
    gens = [mat(_rotation(a, 1)) for a in range(3)]
    I = ((1, 0, 0), (0, 1, 0), (0, 0, 1))

    def produit(A, B):
        return tuple(tuple(sum(A[i][k] * B[k][j] for k in range(3)) for j in range(3)) for i in range(3))

    vus, a_voir = {I}, [I]
    while a_voir:
        x = a_voir.pop()
        for g in gens:
            y = produit(g, x)
            if y not in vus:
                vus.add(y)
                a_voir.append(y)
    return sorted(vus), produit


def appliquer(m, v):
    return tuple(sum(m[i][k] * v[k] for k in range(3)) for i in range(3))


def points_du_cube(genre):
    """Faces (6), sommets (8), arêtes (12) ou grandes diagonales (4) du cube, comme coordonnées."""
    pm = (-1, 1)
    if genre == "faces":
        return [tuple(s if k == a else 0 for k in range(3)) for a in range(3) for s in pm]
    if genre == "sommets":
        return [(x, y, z) for x in pm for y in pm for z in pm]
    if genre == "aretes":
        res = []
        for a in range(3):
            for u in pm:
                for v in pm:
                    autres = iter((u, v))
                    res.append(tuple(0 if k == a else next(autres) for k in range(3)))
        return res
    if genre == "diagonales":
        return [(1, y, z) for y in pm for z in pm]
    raise ValueError(genre)


def action_sur(genre, matrice):
    """Permutation des points d'un genre donné induite par une rotation du cube."""
    pts = points_du_cube(genre)
    idx = {}
    for i, p in enumerate(pts):
        idx[p] = i
        if genre == "diagonales":
            idx[tuple(-x for x in p)] = i          # une diagonale est la même dans les deux sens
    return tuple(idx[appliquer(matrice, p)] for p in pts)
