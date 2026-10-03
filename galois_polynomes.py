"""Le groupe de Galois d'un polynôme à coefficients entiers, calculé.

Deux outils complémentaires :

1. **Degrés 2, 3, 4** : on détermine le groupe exactement avec le discriminant et, pour le
   degré 4, la résolvante cubique (critères classiques, dont celui de Kappe et Warren).

2. **Théorème de Dedekind** : si p ne divise pas le discriminant, la factorisation de f modulo p
   donne le type de cycles d'une permutation du groupe de Galois. Avec beaucoup de nombres
   premiers, on voit apparaître *tous* les types du groupe, dans les proportions prédites par
   le théorème de densité de Chebotarev.

Un polynôme est une liste de coefficients du terme de plus haut degré au terme constant :
x³ − 2 s'écrit [1, 0, 0, -2].
"""

from fractions import Fraction
from math import isqrt

import catalogue as cat
from permgroup import type_de_cycles


# ----------------------------------------------------------------------
# Outils de base
# ----------------------------------------------------------------------
def _nettoyer(f):
    f = list(f)
    while len(f) > 1 and f[0] == 0:
        f.pop(0)
    return f


def evaluer(f, x):
    r = 0
    for c in f:
        r = r * x + c
    return r


def derivee(f):
    n = len(f) - 1
    return [c * (n - i) for i, c in enumerate(f[:-1])] or [0]


def est_carre(n):
    return n >= 0 and isqrt(n) ** 2 == n


def diviseurs(n):
    n = abs(n)
    return [d for d in range(1, n + 1) if n % d == 0]


def determinant(m):
    """Déterminant exact d'une matrice de Fractions (élimination de Gauss)."""
    m = [[Fraction(x) for x in ligne] for ligne in m]
    n, det = len(m), Fraction(1)
    for i in range(n):
        piv = next((r for r in range(i, n) if m[r][i] != 0), None)
        if piv is None:
            return Fraction(0)
        if piv != i:
            m[i], m[piv] = m[piv], m[i]
            det = -det
        det *= m[i][i]
        for r in range(i + 1, n):
            facteur = m[r][i] / m[i][i]
            for c in range(i, n):
                m[r][c] -= facteur * m[i][c]
    return det


def discriminant(f):
    """Discriminant exact (entier) : disc = (−1)^(n(n−1)/2) · Res(f, f') / a_n."""
    f = _nettoyer(f)
    n = len(f) - 1
    if n < 1:
        raise ValueError("degré ≥ 1 requis")
    g = derivee(f)
    taille = 2 * n - 1
    lignes = []
    for i in range(n - 1):
        lignes.append([0] * i + f + [0] * (taille - len(f) - i))
    for i in range(n):
        lignes.append([0] * i + g + [0] * (taille - len(g) - i))
    res = determinant(lignes)
    d = (-1) ** (n * (n - 1) // 2) * res / f[0]
    assert d.denominator == 1
    return int(d)


def racines_rationnelles(f):
    """Racines rationnelles distinctes de f (théorème des racines rationnelles)."""
    f = _nettoyer(f)
    if f[-1] == 0:
        return sorted({Fraction(0)} | set(racines_rationnelles(f[:-1] if len(f) > 1 else [0])))
    cand = {Fraction(s * p, q) for p in diviseurs(f[-1]) for q in diviseurs(f[0]) for s in (1, -1)}
    return sorted(r for r in cand if sum(Fraction(c) * r ** (len(f) - 1 - i) for i, c in enumerate(f)) == 0)


def _deg_quartique_facteur_quadratique(f):
    """Si le quartique unitaire f s'écrit (x²+ax+b)(x²+cx+d) à coefficients entiers, renvoie (a,b,c,d)."""
    _, p, q, r, s = f
    for b in diviseurs(s) + [-d for d in diviseurs(s)] if s else [0]:
        d = s // b if b else 0
        if b * d != s:
            continue
        u = q - b - d                      # a·c = u et a + c = p
        disc = p * p - 4 * u
        if disc < 0 or not est_carre(disc):
            continue
        for a in {(p + isqrt(disc)) // 2, (p - isqrt(disc)) // 2}:
            c = p - a
            if a * c == u and a * d + b * c == r:
                return a, b, c, d
    return None


# ----------------------------------------------------------------------
# Groupe de Galois exact, degrés 2 à 4
# ----------------------------------------------------------------------
def groupe_de_galois(f):
    """Groupe de Galois sur Q d'un polynôme unitaire à coefficients entiers de degré 2, 3 ou 4.

    Retourne un dict : nom, ordre, raisons (phrases justifiant la réponse). Un polynôme de degré
    ≤ 4 est toujours résoluble par radicaux : le groupe est un sous-groupe de S4, qui l'est."""
    f = _nettoyer(f)
    if f[0] != 1:
        raise ValueError("polynôme unitaire attendu (coefficient dominant 1)")
    if not 2 <= len(f) - 1 <= 4:
        raise ValueError("degrés 2, 3 et 4 seulement ; pour le degré 5 voir prouve_s5")
    D = discriminant(f)
    if D == 0:
        raise ValueError("polynôme avec une racine multiple : factorisez-le d'abord")
    raisons = [f"discriminant = {D}" + (" (un carré)" if est_carre(D) else " (pas un carré)")]
    nom, ordre = _galois(f, raisons)
    return {"nom": nom, "ordre": ordre, "raisons": raisons}


def _galois(f, raisons):
    n = len(f) - 1
    D = discriminant(f)
    if n == 1:
        return "trivial", 1
    if n == 2:
        if est_carre(D):
            raisons.append("deux racines rationnelles : rien à permuter")
            return "trivial", 1
        raisons.append("irréductible : le groupe échange les deux racines")
        return "Z/2Z", 2
    entieres = [int(r) for r in racines_rationnelles(f) if r.denominator == 1]
    if entieres:                                    # on divise par (x − r) et on recommence
        r, quotient, acc = entieres[0], [], 0
        for c in f:
            acc = acc * r + c
            quotient.append(acc)
        quotient.pop()
        raisons.append(f"racine rationnelle {r} : on se ramène à un polynôme de degré {n - 1}")
        return _galois(quotient, raisons)
    if n == 3:
        if est_carre(D):
            raisons.append("irréductible, discriminant carré : le groupe est dans A3")
            return "Z/3Z", 3
        raisons.append("irréductible, discriminant non carré : toutes les permutations de S3")
        return "S3", 6
    return _galois_quartique(f, D, raisons)


def _galois_quartique(f, D, raisons):
    fact = _deg_quartique_facteur_quadratique(f)
    if fact:
        a, b, c, d = fact
        d1, d2 = a * a - 4 * b, c * c - 4 * d
        raisons.append(f"se factorise en (x²+{a}x+{b})(x²+{c}x+{d}) de discriminants {d1} et {d2}")
        if est_carre(d1) and est_carre(d2):
            return "trivial", 1
        if est_carre(d1) or est_carre(d2) or est_carre(d1 * d2):
            return "Z/2Z", 2
        return "V4", 4
    _, a, b, c, d = f
    R = [1, -b, a * c - 4 * d, -(a * a * d - 4 * b * d + c * c)]
    racines = [int(r) for r in racines_rationnelles(R) if r.denominator == 1]
    raisons.append(f"irréductible ; résolvante cubique {_polynome(R)}")
    if not racines:
        raisons.append("la résolvante est irréductible")
        return ("A4", 12) if est_carre(D) else ("S4", 24)
    if len(racines) == 3:
        raisons.append("la résolvante a trois racines rationnelles")
        return "V4", 4
    r = racines[0]
    raisons.append(f"la résolvante a exactement une racine rationnelle, {r}")
    d1, d2 = r * r - 4 * d, a * a - 4 * (b - r)
    if _carre_dans_q_racine(d1, D) and _carre_dans_q_racine(d2, D):
        raisons.append("x²−rx+d et x²+ax+b−r se décomposent sur Q(√disc) : cyclique (Kappe–Warren)")
        return "Z/4Z", 4
    raisons.append("ces deux trinômes ne se décomposent pas sur Q(√disc) : diédral (Kappe–Warren)")
    return "D4", 8


def _carre_dans_q_racine(delta, D):
    """delta est-il un carré dans Q(√D) ? ⇔ delta ou delta·D est un carré dans Q."""
    return est_carre(delta) or est_carre(delta * D)


def _polynome(coefs):
    """[1, 0, -3, 2] -> 'x³ − 3x + 2'."""
    n = len(coefs) - 1
    sup = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
    morceaux = []
    for i, c in enumerate(coefs):
        d = n - i
        if c == 0:
            continue
        mono = "" if d == 0 else ("x" if d == 1 else "x" + str(d).translate(sup))
        corps = (str(abs(c)) if (abs(c) != 1 or d == 0) else "") + mono
        morceaux.append(("−" if c < 0 else "+") + " " + corps)
    texte = " ".join(morceaux)
    return texte[2:] if texte.startswith("+ ") else "−" + texte[1:] if texte.startswith("−") else texte


def groupe_catalogue(nom):
    """L'objet ``Groupe`` correspondant à un nom renvoyé par ``groupe_de_galois``."""
    table = {"trivial": lambda: cat.cyclique(1), "Z/2Z": lambda: cat.cyclique(2),
             "Z/3Z": lambda: cat.cyclique(3), "Z/4Z": lambda: cat.cyclique(4),
             "V4": cat.klein, "S3": lambda: cat.symetrique(3), "D4": lambda: cat.diedral(4),
             "A4": lambda: cat.alterne(4), "S4": lambda: cat.symetrique(4)}
    return table[nom]()


# ----------------------------------------------------------------------
# Dedekind : factoriser modulo p
# ----------------------------------------------------------------------
def _norm(f, p):
    f = [c % p for c in f]
    while f and f[0] == 0:
        f.pop(0)
    return f


def _divmod_mod(a, b, p):
    """Division euclidienne dans F_p[x] : (quotient, reste)."""
    a, b = _norm(a, p), _norm(b, p)
    inv = pow(b[0], -1, p)
    q = []
    while len(a) >= len(b):
        c = a[0] * inv % p
        q.append(c)
        for i in range(len(b)):
            a[i] = (a[i] - c * b[i]) % p
        a.pop(0)
    return q or [0], _norm(a, p)


def _mul_mod(a, b, p):
    if not a or not b:
        return []
    res = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            res[i + j] = (res[i + j] + x * y) % p
    return _norm(res, p)


def _soustraire(a, b, p):
    taille = max(len(a), len(b))
    a, b = [0] * (taille - len(a)) + list(a), [0] * (taille - len(b)) + list(b)
    return _norm([(x - y) % p for x, y in zip(a, b)], p)


def _gcd_mod(a, b, p):
    a, b = _norm(a, p), _norm(b, p)
    while b:
        a, b = b, _divmod_mod(a, b, p)[1]
    inv = pow(a[0], -1, p)
    return [c * inv % p for c in a]


def _puissance_x_mod(e, f, p):
    """x^e mod f dans F_p[x], par exponentiation rapide."""
    resultat, base = [1], _divmod_mod([1, 0], f, p)[1]
    while e:
        if e & 1:
            resultat = _divmod_mod(_mul_mod(resultat, base, p), f, p)[1]
        base = _divmod_mod(_mul_mod(base, base, p), f, p)[1]
        e >>= 1
    return resultat


def types_de_cycles_mod_p(f, p):
    """Degrés des facteurs irréductibles de f modulo p, du plus grand au plus petit
    (None si f a une racine multiple modulo p, c'est-à-dire si p divise le discriminant).

    Par le théorème de Dedekind, c'est le type de cycles d'un élément du groupe de Galois.
    Méthode : pgcd(x^(p^k) − x, f) est le produit des facteurs de degré divisant k."""
    f = _norm(f, p)
    df = _norm(derivee(f), p)
    if len(f) < 2 or not df or len(_gcd_mod(f, df, p)) > 1:
        return None
    degres, reste, k = [], f, 1
    while len(reste) - 1 >= 2 * k:
        diff = _soustraire(_puissance_x_mod(p ** k, reste, p), [1, 0], p)
        g = _gcd_mod(reste, diff, p) if diff else reste
        d = len(g) - 1
        if d > 0:
            degres += [k] * (d // k)
            reste = _divmod_mod(reste, g, p)[0]
        k += 1
    if len(reste) - 1 > 0:
        degres.append(len(reste) - 1)
    return tuple(sorted(degres, reverse=True))


def premiers_jusqu_a(n):
    crible = bytearray([1]) * (n + 1)
    crible[:2] = b"\x00\x00"
    for i in range(2, isqrt(n) + 1):
        if crible[i]:
            crible[i * i::i] = bytearray(len(crible[i * i::i]))
    return [i for i, v in enumerate(crible) if v]


def statistique_de_frobenius(f, borne=2000):
    """Pour chaque premier p ≤ borne ne divisant pas le discriminant : le type de cycles mod p.

    Retourne ({type: nombre de premiers}, nombre total de premiers utilisés)."""
    D = discriminant(f)
    compte, total = {}, 0
    for p in premiers_jusqu_a(borne):
        if D % p == 0:
            continue
        t = types_de_cycles_mod_p(f, p)
        if t is None:
            continue
        compte[t] = compte.get(t, 0) + 1
        total += 1
    return dict(sorted(compte.items(), reverse=True)), total


def frequences_theoriques(permutations):
    """Proportion d'éléments de chaque type de cycles dans un groupe (théorème de Chebotarev)."""
    n = len(permutations)
    compte = {}
    for p in permutations:
        t = type_de_cycles(p)
        compte[t] = compte.get(t, 0) + 1
    return {t: c / n for t, c in sorted(compte.items(), reverse=True)}


def prouve_s5(f, borne=2000):
    """Preuve qu'un quintique unitaire entier irréductible a pour groupe de Galois S5.

    Ingrédients (théorème de Dedekind) : un premier où f reste irréductible donne un 5-cycle, donc
    le groupe est transitif d'ordre divisible par 5 ; un premier où f a un facteur quadratique et
    trois facteurs linéaires donne une transposition. Un sous-groupe transitif de S_p (p premier)
    contenant une transposition est S_p. Retourne (booléen, p_5cycle, p_transposition)."""
    D = discriminant(f)
    p5 = ptr = None
    for p in premiers_jusqu_a(borne):
        if D % p == 0:
            continue
        t = types_de_cycles_mod_p(f, p)
        if t == (5,) and p5 is None:
            p5 = p
        if t == (2, 1, 1, 1) and ptr is None:
            ptr = p
        if p5 and ptr:
            return True, p5, ptr
    return False, p5, ptr
