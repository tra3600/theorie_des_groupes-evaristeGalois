"""Noyau mathématique : la classe ``Groupe`` et ses outils d'exploration.

Un groupe fini est décrit par la liste de ses éléments et une opération
binaire.  À la construction, on calcule une fois pour toutes la table de
Cayley (sous forme d'indices entiers) : toutes les vérifications ensuite
(associativité, sous-groupes, série dérivée...) sont de simples lectures
dans cette table, ce qui permet de manipuler S5 (120 éléments) sans attendre.
"""

from dataclasses import dataclass, field
from itertools import product


@dataclass
class Axiome:
    """Résultat de la vérification d'un axiome de groupe."""
    nom: str
    valide: bool
    explication: str


@dataclass
class Rapport:
    """Bilan complet : le couple (ensemble, opération) est-il un groupe ?"""
    axiomes: list = field(default_factory=list)

    @property
    def est_groupe(self):
        return all(a.valide for a in self.axiomes)

    def __str__(self):
        lignes = []
        for a in self.axiomes:
            marque = "✔" if a.valide else "✘"
            lignes.append(f"  {marque} {a.nom:<14} {a.explication}")
        verdict = "C'est un groupe !" if self.est_groupe else "Ce n'est PAS un groupe."
        lignes.append(f"  → {verdict}")
        return "\n".join(lignes)


class Groupe:
    """Un groupe fini (ou un candidat au titre de groupe).

    Paramètres
    ----------
    elements : itérable d'objets hashables
    operation : fonction (a, b) -> a·b
    nom : nom affiché du groupe
    etiquette : fonction élément -> texte lisible (par défaut ``str``)
    """

    def __init__(self, elements, operation, nom="G", etiquette=str):
        self.elements = list(elements)
        self.operation = operation
        self.nom = nom
        self._etiquette = etiquette
        self._index = {e: i for i, e in enumerate(self.elements)}
        if len(self._index) != len(self.elements):
            raise ValueError("Les éléments d'un groupe doivent être distincts.")
        # table[i][j] = indice de elements[i]·elements[j], ou None si le
        # résultat sort de l'ensemble (fermeture violée).
        self.table = [
            [self._index.get(operation(a, b)) for b in self.elements]
            for a in self.elements
        ]
        self._neutre = self._chercher_neutre()
        self._inv = self._chercher_inverses()

    # ------------------------------------------------------------------
    # Accès de base
    # ------------------------------------------------------------------
    def __len__(self):
        return len(self.elements)

    def __iter__(self):
        return iter(self.elements)

    def __contains__(self, x):
        return x in self._index

    def __repr__(self):
        return f"<Groupe {self.nom} d'ordre {len(self)}>"

    def etiquette(self, e):
        return self._etiquette(e)

    def op(self, a, b):
        """Produit a·b (lève une erreur si le résultat sort de l'ensemble)."""
        k = self.table[self._index[a]][self._index[b]]
        if k is None:
            raise ValueError(f"{self.etiquette(a)}·{self.etiquette(b)} sort de l'ensemble")
        return self.elements[k]

    def produit(self, *xs):
        """Produit de plusieurs éléments, de gauche à droite."""
        r = self.neutre
        for x in xs:
            r = self.op(r, x)
        return r

    @property
    def neutre(self):
        return None if self._neutre is None else self.elements[self._neutre]

    def inverse(self, a):
        k = self._inv.get(self._index[a])
        return None if k is None else self.elements[k]

    # Compatibilité avec l'ancienne API (attributs en anglais)
    @property
    def identity(self):
        return self.neutre

    @property
    def inverses(self):
        return {self.elements[i]: self.elements[j] for i, j in self._inv.items()}

    def is_group(self):
        return self.verifier().est_groupe

    # ------------------------------------------------------------------
    # Axiomes
    # ------------------------------------------------------------------
    def _chercher_neutre(self):
        n = len(self)
        for e in range(n):
            if all(self.table[e][a] == a and self.table[a][e] == a for a in range(n)):
                return e
        return None

    def _chercher_inverses(self):
        inv = {}
        e = self._neutre
        if e is None:
            return inv
        for a in range(len(self)):
            for b in range(len(self)):
                if self.table[a][b] == e and self.table[b][a] == e:
                    inv[a] = b
                    break
        return inv

    def verifier(self):
        """Vérifie les 4 axiomes et explique, avec un contre-exemple, ce qui cloche."""
        E, T, lab = self.elements, self.table, self.etiquette
        n = len(self)
        rapport = Rapport()

        # 1. Fermeture (loi interne)
        fermeture = Axiome("Fermeture", True, "a·b reste toujours dans l'ensemble")
        for i, j in product(range(n), repeat=2):
            if T[i][j] is None:
                resultat = self.operation(E[i], E[j])
                fermeture = Axiome(
                    "Fermeture", False,
                    f"{lab(E[i])}·{lab(E[j])} = {resultat} n'est pas dans l'ensemble")
                break
        rapport.axiomes.append(fermeture)

        # 2. Associativité (seulement testable si la loi est interne)
        if not fermeture.valide:
            rapport.axiomes.append(Axiome("Associativité", False, "non testable sans fermeture"))
        else:
            assoc = Axiome("Associativité", True, "(a·b)·c = a·(b·c) pour tous a, b, c")
            for i, j in product(range(n), repeat=2):
                ij = T[i][j]
                ligne_ij, ligne_i = T[ij], T[i]
                for k in range(n):
                    if ligne_ij[k] != ligne_i[T[j][k]]:
                        assoc = Axiome(
                            "Associativité", False,
                            f"({lab(E[i])}·{lab(E[j])})·{lab(E[k])} = {lab(E[ligne_ij[k]])}"
                            f" mais {lab(E[i])}·({lab(E[j])}·{lab(E[k])}) = "
                            f"{lab(E[ligne_i[T[j][k]]])}")
                        break
                if not assoc.valide:
                    break
            rapport.axiomes.append(assoc)

        # 3. Élément neutre
        if self._neutre is None:
            rapport.axiomes.append(Axiome("Neutre", False, "aucun élément e tel que e·a = a·e = a"))
        else:
            rapport.axiomes.append(Axiome("Neutre", True, f"e = {lab(self.neutre)}"))

        # 4. Inverses
        if self._neutre is None:
            rapport.axiomes.append(Axiome("Inverses", False, "impossible sans élément neutre"))
        else:
            sans_inverse = [E[a] for a in range(n) if a not in self._inv]
            if sans_inverse:
                rapport.axiomes.append(Axiome(
                    "Inverses", False,
                    f"{lab(sans_inverse[0])} n'a pas d'inverse"
                    + (f", ni {len(sans_inverse) - 1} autre(s) élément(s)" if len(sans_inverse) > 1 else "")))
            else:
                rapport.axiomes.append(Axiome("Inverses", True, "chaque élément a un symétrique"))
        return rapport

    # ------------------------------------------------------------------
    # Éléments
    # ------------------------------------------------------------------
    def puissance(self, a, k):
        """a^k (k peut être négatif)."""
        if k < 0:
            a, k = self.inverse(a), -k
        r = self.neutre
        for _ in range(k):
            r = self.op(r, a)
        return r

    def ordre_element(self, a):
        """Plus petit k ≥ 1 tel que a^k = e."""
        e, i, x, k = self._neutre, self._index[a], self._index[a], 1
        while x != e:
            x = self.table[x][i]
            k += 1
        return k

    def est_abelien(self):
        return self.contre_exemple_commutativite() is None

    def contre_exemple_commutativite(self):
        """Renvoie (a, b) avec a·b ≠ b·a, ou None si le groupe est commutatif."""
        n = len(self)
        for i in range(n):
            for j in range(i + 1, n):
                if self.table[i][j] != self.table[j][i]:
                    return self.elements[i], self.elements[j]
        return None

    def generateur(self):
        """Un élément qui engendre tout le groupe, ou None (groupe non cyclique)."""
        for a in self.elements:
            if self.ordre_element(a) == len(self):
                return a
        return None

    def est_cyclique(self):
        return self.generateur() is not None

    # ------------------------------------------------------------------
    # Sous-groupes
    # ------------------------------------------------------------------
    def _cloture(self, indices):
        """Plus petit sous-groupe (ensemble d'indices) contenant ``indices``."""
        sg = {self._neutre} | set(indices)
        a_traiter = list(sg)
        while a_traiter:
            x = a_traiter.pop()
            for y in list(sg):
                for z in (self.table[x][y], self.table[y][x]):
                    if z not in sg:
                        sg.add(z)
                        a_traiter.append(z)
        return frozenset(sg)

    def engendre(self, *generateurs):
        """Sous-groupe engendré, renvoyé comme liste d'éléments."""
        idx = self._cloture(self._index[g] for g in generateurs)
        return [self.elements[i] for i in sorted(idx)]

    def sous_groupe(self, elements, nom="H"):
        """Construit un objet ``Groupe`` à partir d'une partie stable."""
        return Groupe(elements, self.operation, nom=nom, etiquette=self._etiquette)

    def sous_groupes(self):
        """Tous les sous-groupes (listes d'éléments), triés par ordre croissant.

        Méthode : on part des sous-groupes cycliques puis on combine deux à
        deux jusqu'à stabilisation.  Parfait pour les groupes de taille < 50.
        """
        cycliques = {self._cloture([i]) for i in range(len(self))}
        tous = set(cycliques)
        frontiere = set(cycliques)
        while frontiere:
            nouveaux = set()
            for H in frontiere:
                for C in cycliques:
                    if not C <= H:
                        K = self._cloture(H | C)
                        if K not in tous:
                            nouveaux.add(K)
            tous |= nouveaux
            frontiere = nouveaux
        return [[self.elements[i] for i in sorted(H)]
                for H in sorted(tous, key=lambda H: (len(H), sorted(H)))]

    def classes_a_gauche(self, H):
        """Les classes aH : elles découpent le groupe en morceaux de même taille."""
        h = [self._index[x] for x in H]
        vues, classes = set(), []
        for a in range(len(self)):
            if a not in vues:
                classe = sorted({self.table[a][x] for x in h})
                vues.update(classe)
                classes.append([self.elements[i] for i in classe])
        return classes

    def est_distingue(self, H):
        """H est distingué (normal) si gHg⁻¹ = H pour tout g."""
        h = {self._index[x] for x in H}
        for g in range(len(self)):
            gi = self._inv[g]
            if any(self.table[self.table[g][x]][gi] not in h for x in h):
                return False
        return True

    def conjugue(self, g, x):
        """g·x·g⁻¹"""
        return self.produit(g, x, self.inverse(g))

    def classes_de_conjugaison(self):
        vues, classes = set(), []
        for x in self.elements:
            if x not in vues:
                classe = []
                for g in self.elements:
                    y = self.conjugue(g, x)
                    if y not in classe:
                        classe.append(y)
                vues.update(classe)
                classes.append(classe)
        return classes

    # ------------------------------------------------------------------
    # Centre, quotients, produits, homomorphismes
    # ------------------------------------------------------------------
    def centre(self):
        """Z(G) : les éléments qui commutent avec tous les autres."""
        n = len(self)
        return [self.elements[i] for i in range(n)
                if all(self.table[i][j] == self.table[j][i] for j in range(n))]

    def centralisateur(self, a):
        """C(a) : les éléments qui commutent avec a."""
        i = self._index[a]
        return [self.elements[j] for j in range(len(self))
                if self.table[i][j] == self.table[j][i]]

    def equation_des_classes(self):
        """|G| = |Z(G)| + Σ [G : C(a)] sur un représentant a de chaque classe non centrale.

        Retourne (taille du centre, [tailles des classes non triviales])."""
        classes = self.classes_de_conjugaison()
        return len(self.centre()), sorted(len(c) for c in classes if len(c) > 1)

    def normalisateur(self, H):
        """N(H) = {g : gHg⁻¹ = H} (le plus grand sous-groupe où H est distingué)."""
        h = {self._index[x] for x in H}
        res = []
        for g in range(len(self)):
            gi = self._inv[g]
            if {self.table[self.table[g][x]][gi] for x in h} == h:
                res.append(self.elements[g])
        return res

    def quotient(self, N, nom=None):
        """G/N : le groupe des classes aN, pour un sous-groupe distingué N."""
        if not self.est_distingue(N):
            raise ValueError("Le quotient n'existe que pour un sous-groupe distingué.")
        classes = [frozenset(c) for c in self.classes_a_gauche(N)]
        trouver = {x: c for c in classes for x in c}

        def produit(c1, c2):
            return trouver[self.op(next(iter(c1)), next(iter(c2)))]

        def etiquette(c):
            return "{" + ",".join(self.etiquette(x) for x in sorted(c, key=self._index.get)) + "}"

        return Groupe(classes, produit, nom=nom or f"{self.nom}/N", etiquette=etiquette)

    def produit_direct(self, autre, nom=None):
        """G × H, avec l'opération composante par composante."""
        elems = [(a, b) for a in self.elements for b in autre.elements]
        return Groupe(elems, lambda x, y: (self.op(x[0], y[0]), autre.op(x[1], y[1])),
                      nom=nom or f"{self.nom}×{autre.nom}",
                      etiquette=lambda x: f"({self.etiquette(x[0])},{autre.etiquette(x[1])})")

    def est_homomorphisme(self, f, autre):
        """f (dict ou fonction) vérifie-t-elle f(a·b) = f(a)·f(b) vers ``autre`` ?"""
        g = f if callable(f) else f.__getitem__
        return all(g(self.op(a, b)) == autre.op(g(a), g(b))
                   for a in self.elements for b in self.elements)

    def noyau(self, f, autre):
        """{a : f(a) = neutre de l'arrivée}."""
        g = f if callable(f) else f.__getitem__
        return [a for a in self.elements if g(a) == autre.neutre]

    def image(self, f):
        g = f if callable(f) else f.__getitem__
        vus = []
        for a in self.elements:
            y = g(a)
            if y not in vus:
                vus.append(y)
        return vus

    def est_simple(self):
        """Simple : seuls {e} et G sont distingués (on regarde la clôture normale de chaque classe)."""
        if len(self) == 1:
            return False
        for classe in self.classes_de_conjugaison():
            if classe == [self.neutre]:
                continue
            if len(self._cloture(self._index[x] for x in classe)) != len(self):
                return False
        return True

    def sylow(self, p):
        """Un p-sous-groupe de Sylow P et le nombre n_p de ses conjugués.

        P est construit pas à pas : tant que |P| < p^a, on cherche dans le normalisateur de P
        un élément g hors de P tel que <P, g> soit encore un p-groupe de taille p·|P|."""
        n, a = len(self), 0
        while n % p == 0:
            n //= p
            a += 1
        cible = p ** a
        P = [self.neutre]
        while len(P) < cible:
            for g in self.normalisateur(P):
                if g in P:
                    continue
                Q = self.engendre(*P, g)
                if len(Q) == p * len(P):
                    P = Q
                    break
            else:
                raise RuntimeError("construction de Sylow échouée")
        conjugues = {frozenset(self.conjugue(g, x) for x in P) for g in self.elements}
        return P, len(conjugues)

    def representation_reguliere(self):
        """Théorème de Cayley : chaque g agit sur G par x ↦ g·x, c'est une permutation des indices.

        Retourne {g: permutation (tuple d'indices)} ; g ↦ perm(g) est un homomorphisme injectif G → S_|G|."""
        return {g: tuple(self.table[self._index[g]][j] for j in range(len(self))) for g in self.elements}

    # ------------------------------------------------------------------
    # Résolubilité (le cœur de la théorie de Galois)
    # ------------------------------------------------------------------
    def commutateur(self, a, b):
        """[a, b] = a·b·a⁻¹·b⁻¹ : il vaut e exactement quand a et b commutent."""
        return self.produit(a, b, self.inverse(a), self.inverse(b))

    def derive(self, H=None):
        """Sous-groupe dérivé [H, H] engendré par les commutateurs de H."""
        H = self.elements if H is None else H
        comms = {self._index[self.commutateur(a, b)] for a in H for b in H}
        return [self.elements[i] for i in sorted(self._cloture(comms))]

    def serie_derivee(self):
        """G ⊇ G' ⊇ G'' ⊇ ... jusqu'à stabilisation."""
        serie = [list(self.elements)]
        while True:
            suivant = self.derive(serie[-1])
            if len(suivant) == len(serie[-1]):
                return serie
            serie.append(suivant)

    def est_resoluble(self):
        """Résoluble ⇔ la série dérivée finit sur {e}."""
        return len(self.serie_derivee()[-1]) == 1

    # ------------------------------------------------------------------
    # Isomorphisme (recherche simple, pour les petits groupes)
    # ------------------------------------------------------------------
    def signature(self):
        """Nombre d'éléments de chaque ordre : deux groupes isomorphes ont la même."""
        compte = {}
        for a in self.elements:
            k = self.ordre_element(a)
            compte[k] = compte.get(k, 0) + 1
        return dict(sorted(compte.items()))

    def isomorphisme(self, autre):
        """Cherche une bijection f avec f(a·b) = f(a)·f(b), ou renvoie None."""
        if len(self) != len(autre) or self.signature() != autre.signature():
            return None
        n = len(self)
        # On n'a besoin d'envoyer que des générateurs : on en choisit peu.
        gens = []
        for i in sorted(range(n), key=lambda i: -self.ordre_element(self.elements[i])):
            if len(self._cloture(gens)) == n:
                break
            if i not in self._cloture(gens):
                gens.append(i)
        ordres = [self.ordre_element(self.elements[g]) for g in gens]
        candidats = [[j for j in range(n) if autre.ordre_element(autre.elements[j]) == o]
                     for o in ordres]

        def essayer(images):
            f = {self._neutre: autre._neutre}
            f.update(zip(gens, images))
            a_traiter = list(f)
            while a_traiter:
                x = a_traiter.pop()
                for g, fg in zip(gens, images):
                    y, fy = self.table[x][g], autre.table[f[x]][fg]
                    if y in f:
                        if f[y] != fy:
                            return None
                    else:
                        f[y] = fy
                        a_traiter.append(y)
            if len(set(f.values())) != n:
                return None
            if all(f[self.table[a][b]] == autre.table[f[a]][f[b]]
                   for a in range(n) for b in range(n)):
                return {self.elements[a]: autre.elements[b] for a, b in f.items()}
            return None

        for images in product(*candidats):
            f = essayer(images)
            if f:
                return f
        return None

    # ------------------------------------------------------------------
    # Affichage texte
    # ------------------------------------------------------------------
    def table_texte(self, max_elements=12):
        """Table de Cayley en ASCII (tronquée si le groupe est trop gros)."""
        if len(self) > max_elements:
            return f"(table {len(self)}×{len(self)} trop grande pour l'écran)"
        labs = [self.etiquette(e) for e in self.elements]
        w = max(len(s) for s in labs) + 1
        tete = " " * w + "│" + "".join(s.rjust(w) for s in labs)
        lignes = [tete, "─" * w + "┼" + "─" * (w * len(labs))]
        for i, s in enumerate(labs):
            cellules = "".join(
                (labs[k] if k is not None else "?").rjust(w) for k in self.table[i])
            lignes.append(s.rjust(w) + "│" + cellules)
        return "\n".join(lignes)
