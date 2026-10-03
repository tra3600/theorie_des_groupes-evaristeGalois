"""Les groupes (Z/nZ)* au travail : Fermat, Euler, Diffie–Hellman, RSA, logarithme discret.

Tout repose sur le théorème de Lagrange : dans un groupe d'ordre m, a^m = e.
Pour (Z/nZ)*, d'ordre φ(n), cela donne a^φ(n) ≡ 1 (mod n) : le théorème d'Euler.
"""

import random
from math import gcd, isqrt


def phi(n):
    """Indicatrice d'Euler : le nombre d'entiers de 1 à n premiers avec n."""
    r, m, p = n, n, 2
    while p * p <= m:
        if m % p == 0:
            while m % p == 0:
                m //= p
            r -= r // p
        p += 1
    if m > 1:
        r -= r // m
    return r


def facteurs_premiers(n):
    f, p = {}, 2
    while p * p <= n:
        while n % p == 0:
            f[p] = f.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def est_premier(n):
    return n > 1 and facteurs_premiers(n) == {n: 1}


def ordre_multiplicatif(a, n):
    """Ordre de a dans (Z/nZ)* (a premier avec n)."""
    if gcd(a, n) != 1:
        raise ValueError("a doit être premier avec n")
    k, x = 1, a % n
    while x != 1 % n:
        x = x * a % n
        k += 1
    return k


def racines_primitives(p):
    """Les générateurs de (Z/pZ)* (p premier) : les éléments d'ordre p − 1."""
    if not est_premier(p):
        raise ValueError("p doit être premier")
    return [g for g in range(1, p) if ordre_multiplicatif(g, p) == p - 1]


def a_une_racine_primitive(n):
    """(Z/nZ)* est cyclique exactement pour n = 1, 2, 4, p^k, 2p^k (p premier impair)."""
    if n in (1, 2, 4):
        return True
    m = n // 2 if n % 2 == 0 else n
    if n % 2 == 0 and m % 2 == 0:
        return False
    f = facteurs_premiers(m)
    return len(f) == 1 and 2 not in f


def carmichael(n):
    """λ(n) : le plus petit k tel que a^k ≡ 1 (mod n) pour tout a premier avec n (l'exposant du groupe)."""
    from math import lcm
    r = 1
    for p, e in facteurs_premiers(n).items():
        l = (p - 1) * p ** (e - 1) if p > 2 else (1 if e == 1 else 2 if e == 2 else 2 ** (e - 2))
        r = lcm(r, l)
    return r


def echange_diffie_hellman(p, g, secret_a, secret_b):
    """Alice envoie A = g^a, Bob envoie B = g^b ; tous deux obtiennent g^(ab) sans l'avoir échangé.

    Retourne (A, B, clé d'Alice, clé de Bob)."""
    A, B = pow(g, secret_a, p), pow(g, secret_b, p)
    return A, B, pow(B, secret_a, p), pow(A, secret_b, p)


def logarithme_discret(g, h, p):
    """x tel que g^x ≡ h (mod p) par « pas de bébé, pas de géant » : O(√p) au lieu de O(p).

    Retourne (x, nombre d'opérations) ou (None, opérations) si h n'est pas une puissance de g."""
    ordre = ordre_multiplicatif(g, p)
    m = isqrt(ordre) + 1
    bebes = {}
    x, ops = 1, 0
    for j in range(m):
        bebes.setdefault(x, j)
        x = x * g % p
        ops += 1
    facteur = pow(g, -m, p)
    gamma = h % p
    for i in range(m):
        if gamma in bebes:
            return i * m + bebes[gamma], ops
        gamma = gamma * facteur % p
        ops += 1
    return None, ops


def logarithme_naif(g, h, p):
    """Essaie g^0, g^1, ... : (x, nombre d'essais)."""
    x, k = 1, 0
    while x != h % p:
        x = x * g % p
        k += 1
        if k > p:
            return None, k
    return k, k + 1


def rsa_jouet(p, q, e, message):
    """Chiffre puis déchiffre un message avec RSA. Retourne (n, φ(n), d, chiffré, déchiffré).

    d est l'inverse de e dans (Z/φ(n)Z)* ; le déchiffrement marche car, dans (Z/nZ)* d'ordre φ(n),
    m^(e·d) = m^(1 + kφ) = m."""
    n, f = p * q, (p - 1) * (q - 1)
    if gcd(e, f) != 1:
        raise ValueError("e doit être premier avec φ(n)")
    d = pow(e, -1, f)
    c = pow(message, e, n)
    return n, f, d, c, pow(c, d, n)


def nombre_premier_aleatoire(bits, rng=None):
    rng = rng or random.Random(0)
    while True:
        n = rng.getrandbits(bits) | (1 << (bits - 1)) | 1
        if est_premier(n):
            return n


def plus_petite_racine_primitive(p):
    """Le plus petit générateur de (Z/pZ)* (p premier), sans énumérer le groupe :
    g est un générateur ssi g^((p−1)/q) ≠ 1 pour tout facteur premier q de p − 1."""
    qs = list(facteurs_premiers(p - 1))
    return next(g for g in range(2, p) if all(pow(g, (p - 1) // q, p) != 1 for q in qs))
