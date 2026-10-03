"""Le Défi de Galois : un quiz à questions tirées au hasard.

Chaque question est calculée par la bibliothèque elle-même : la réponse
attendue n'est jamais écrite à la main, elle est *démontrée* par le programme.
"""

import random
from dataclasses import dataclass

import arithmetique as ar
import denombrement as dn
import galois_polynomes as gp
import permgroup as pg
from catalogue import (alterne, cyclique, depuis_cycles, diedral, horloge, inversibles_mod,
                       klein, melange_parfait, groupe_engendre_par_permutation,
                       notation_cyclique, quaternions, symetrique)


@dataclass
class Question:
    enonce: str
    reponse: str
    explication: str
    indice: str = ""


def _oui_non(b):
    return "oui" if b else "non"


def _normaliser(texte):
    t = texte.strip().lower().replace(" ", "").replace("h", "")
    return {"o": "oui", "y": "oui", "yes": "oui", "n": "non", "no": "non"}.get(t, t)


# ----------------------------------------------------------------------
# Générateurs de questions
# ----------------------------------------------------------------------
def q_horloge(rng):
    G = horloge()
    a, b = rng.randint(1, 12), rng.randint(5, 30)
    r = G.op(a % 12, b % 12)
    return Question(
        f"Il est {a}h. Quelle heure affichera l'horloge dans {b} heures ?",
        str(r or 12),
        f"{a} + {b} = {a + b} ≡ {r or 12} (mod 12) : on ne garde que le reste de la division par 12.",
        "Enlève des paquets de 12 heures.")


def q_inverse_modulo(rng):
    p = rng.choice([5, 7, 11, 13])
    G = inversibles_mod(p)
    a = rng.randint(2, p - 1)
    inv = G.inverse(a)
    return Question(
        f"Dans (Z/{p}Z)*, quel nombre x vérifie {a} × x ≡ 1 (mod {p}) ?",
        str(inv),
        f"{a} × {inv} = {a * inv} = {a * inv // p}×{p} + 1.",
        f"Essaie x = 2, 3, ... jusqu'à {p - 1}.")


def q_ordre_permutation(rng):
    n = rng.choice([4, 5, 6, 7])
    p = list(range(n))
    rng.shuffle(p)
    p = tuple(p)
    if p == tuple(range(n)):
        p = depuis_cycles(n, (1, 2))
    G = groupe_engendre_par_permutation(p)
    k = len(G)
    return Question(
        f"Combien de fois faut-il répéter la permutation {notation_cyclique(p)} "
        f"pour revenir à la position de départ ?",
        str(k),
        "L'ordre d'une permutation est le PPCM des longueurs de ses cycles.",
        "Regarde la longueur de chaque cycle.")


def q_commutatif(rng):
    G = rng.choice([cyclique(6), klein(), symetrique(3), diedral(4), quaternions(),
                    inversibles_mod(8), diedral(5), cyclique(8)])
    ce = G.contre_exemple_commutativite()
    if ce:
        a, b = ce
        expl = (f"Non : {G.etiquette(a)}·{G.etiquette(b)} = {G.etiquette(G.op(a, b))} mais "
                f"{G.etiquette(b)}·{G.etiquette(a)} = {G.etiquette(G.op(b, a))}.")
    else:
        expl = "Oui : dans sa table de Cayley, la case (a, b) égale toujours la case (b, a)."
    return Question(f"Le groupe {G.nom} ({len(G)} éléments) est-il commutatif ? (oui/non)",
                    _oui_non(ce is None), expl,
                    "Cherche deux éléments a, b avec a·b ≠ b·a (pense aux polygones et aux permutations).")


def q_cardinal_sn(rng):
    n = rng.choice([3, 4, 5, 6])
    import math
    return Question(
        f"De combien de façons peut-on ranger {n} objets différents ? (cardinal de S{n})",
        str(math.factorial(n)),
        f"{n}! = " + " × ".join(str(i) for i in range(n, 0, -1)) + f" = {math.factorial(n)}.",
        "n choix pour le premier, n−1 pour le deuxième...")


def q_symetries_polygone(rng):
    n = rng.randint(3, 9)
    return Question(
        f"Combien de symétries possède un polygone régulier à {n} côtés ?",
        str(len(diedral(n))),
        f"{n} rotations + {n} miroirs = {2 * n} : c'est le groupe diédral D{n}.",
        "Compte les rotations, puis les axes de symétrie.")


def q_composition_polygone(rng):
    n = rng.choice([3, 4])
    G = diedral(n)
    a, b = rng.sample(G.elements, 2)
    return Question(
        f"Dans D{n} (symétries du {'triangle' if n == 3 else 'carré'}), "
        f"que vaut {G.etiquette(a)} ∘ {G.etiquette(b)} ? (réponse du type R90°, M2, id)",
        G.etiquette(G.op(a, b)),
        "On applique d'abord la symétrie de droite, puis celle de gauche. "
        "Deux miroirs donnent une rotation ; un miroir et une rotation donnent un miroir.",
        "Lance l'illustration « triangle » ou « carré » pour voir les dessins.")


def q_melange(rng):
    n = rng.choice([4, 6, 8, 10, 52])
    k = len(groupe_engendre_par_permutation(melange_parfait(n)))
    return Question(
        f"Avec un paquet de {n} cartes, combien de mélanges parfaits (faro extérieur) "
        f"faut-il pour retrouver l'ordre initial ?",
        str(k),
        f"Le mélange est une permutation d'ordre {k} dans S{n}.",
        "Pour 52 cartes, les magiciens connaissent la réponse par cœur !")


def q_lagrange(rng):
    G = rng.choice([symetrique(3), diedral(4), cyclique(12), quaternions(), symetrique(4)])
    diviseurs = [d for d in range(1, len(G) + 1) if len(G) % d == 0]
    intrus = rng.choice([d for d in range(2, len(G)) if len(G) % d] or [len(G) + 1])
    candidats = sorted(rng.sample(diviseurs, min(2, len(diviseurs))) + [intrus])
    return Question(
        f"{G.nom} a {len(G)} éléments. Lequel de ces nombres NE PEUT PAS être la taille "
        f"d'un de ses sous-groupes : {', '.join(map(str, candidats))} ?",
        str(intrus),
        f"Théorème de Lagrange : la taille d'un sous-groupe divise celle du groupe, "
        f"et {intrus} ne divise pas {len(G)}.",
        "Pense à la divisibilité.")


def q_sylow(rng):
    G = rng.choice([symetrique(4), alterne(4), diedral(6), alterne(5)])
    p = rng.choice(list(ar.facteurs_premiers(len(G))))
    P, n = G.sylow(p)
    return Question(
        f"Combien {G.nom} a-t-il de sous-groupes de Sylow {p} (de taille {len(P)}) ?",
        str(n),
        f"Le programme en construit un et compte ses conjugués : {n}. "
        f"Sylow : {n} ≡ 1 (mod {p}) et {n} divise {len(G) // len(P)}.",
        f"n_{p} est congru à 1 modulo {p} et divise {len(G) // len(P)}.")


def q_burnside(rng):
    if rng.random() < 0.5:
        n, k = rng.choice([4, 5, 6, 7]), rng.choice([2, 3])
        r = dn.colliers(n, k)
        return Question(
            f"Combien de colliers différents avec {n} perles de {k} couleurs (rotation seulement) ?",
            str(r),
            f"Burnside : moyenne de k^(cycles) sur les {n} rotations = {r}.",
            f"Il y a {k}^{n} = {k ** n} coloriages numérotés ; un collier en regroupe plusieurs.")
    k = rng.choice([2, 3])
    r = dn.coloriages_du_cube("faces", k)
    return Question(
        f"Combien de cubes différents peut-on peindre en colorant chaque face de l'une des {k} couleurs "
        f"(à rotation près) ?", str(r),
        f"Burnside sur les 24 rotations du cube : {r}.", "Il y a 24 rotations ; compte leurs cycles sur les 6 faces.")


def q_diffie_hellman(rng):
    p = rng.choice([11, 13, 17, 19, 23])
    g = rng.choice(ar.racines_primitives(p))
    a, b = rng.randint(2, p - 2), rng.randint(2, p - 2)
    A, B, cle, _ = ar.echange_diffie_hellman(p, g, a, b)
    return Question(
        f"Diffie–Hellman avec p = {p}, g = {g} : Alice envoie {A} = g^a et son secret est a = {a}. "
        f"Bob envoie {B}. Quelle est la clé commune ?", str(cle),
        f"{B}^{a} mod {p} = {cle}  (et {A}^b = g^(ab) donne la même chose pour Bob).",
        f"Calcule {B}^{a} modulo {p}.")


def q_rsa(rng):
    p, q = rng.choice([(5, 11), (7, 11), (7, 13), (11, 13), (5, 17)])
    f = (p - 1) * (q - 1)
    return Question(
        f"RSA avec p = {p} et q = {q} : combien vaut φ(n), l'ordre du groupe (Z/nZ)* ?",
        str(f), f"φ({p * q}) = ({p} − 1)({q} − 1) = {f}.", "Pour n = p·q, φ(n) = (p − 1)(q − 1).")


def q_galois_polynome(rng):
    nom, f = rng.choice([("x³ − 2", [1, 0, 0, -2]), ("x³ − 3x + 1", [1, 0, -3, 1]),
                         ("x³ − x − 1", [1, 0, -1, -1]), ("x⁴ − 2", [1, 0, 0, 0, -2]),
                         ("x⁴ + 1", [1, 0, 0, 0, 1]), ("x⁴ − x − 1", [1, 0, 0, -1, -1]),
                         ("x⁴ + x³ + x² + x + 1", [1, 1, 1, 1, 1])])
    r = gp.groupe_de_galois(f)
    return Question(
        f"Quel est l'ordre du groupe de Galois de {nom} sur Q ?", str(r["ordre"]),
        f"C'est {r['nom']}. " + " ; ".join(r["raisons"]) + ".",
        "Regarde le discriminant (carré ou non ?) et la résolvante cubique au degré 4.")


def q_centre(rng):
    n = rng.choice([3, 4, 5, 6, 7, 8])
    G = diedral(n)
    return Question(
        f"Combien d'éléments le centre de D{n} (les 2·{n} symétries du {n}-gone) contient-il ?",
        str(len(G.centre())),
        f"Le centre est {{{', '.join(G.etiquette(x) for x in G.centre())}}} : "
        + ("seule la rotation d'un demi-tour commute avec tout quand n est pair." if n % 2 == 0
           else "pour n impair, seule l'identité commute avec tout."),
        "Quelles symétries commutent avec chaque miroir ?")


def q_rubik(rng):
    R, U = pg.mouvement_cube("R"), pg.mouvement_cube("U")
    Ri, Ui = pg.inverse(R), pg.inverse(U)
    nom, p = rng.choice([("R", R), ("R U", pg.compose(R, U)),
                         ("R U R' U'", pg.compose(pg.compose(R, U), pg.compose(Ri, Ui)))])
    o = pg.ordre_permutation(p)
    return Question(
        f"Sur un Rubik's Cube, combien de fois faut-il répéter l'algorithme « {nom} » pour retrouver "
        f"le cube de départ ?", str(o),
        f"L'ordre de cette permutation des autocollants est {o} (ppcm des longueurs de ses cycles).",
        "Calcule le ppcm des longueurs des cycles.")


def q_simple(rng):
    choix = [("Z/6Z", cyclique(6)), ("S4", symetrique(4)), ("A4", alterne(4)), ("D5", diedral(5)),
             ("A5", alterne(5)), ("Z/7Z", cyclique(7))]
    rng.shuffle(choix)
    proposes = choix[:4]
    simples = [n for n, G in proposes if G.est_simple()]
    if len(simples) != 1:
        proposes = [c for c in choix if c[0] == "A5"] + [c for c in choix if c[0] in ("S4", "A4", "D5")][:3]
        simples = ["A5"]
    return Question(
        f"Lequel de ces groupes est simple (sans sous-groupe distingué autre que {{e}} et lui-même) : "
        f"{', '.join(n for n, _ in proposes)} ?", simples[0],
        f"{simples[0]} est le seul simple ; les autres ont un sous-groupe distingué non trivial.",
        "Un groupe cyclique d'ordre non premier a des sous-groupes ; A5 est le premier simple non abélien.")


QUESTIONS = [q_horloge, q_inverse_modulo, q_ordre_permutation, q_commutatif,
             q_cardinal_sn, q_symetries_polygone, q_composition_polygone,
             q_melange, q_lagrange, q_sylow, q_burnside, q_diffie_hellman, q_rsa,
             q_galois_polynome, q_centre, q_rubik, q_simple]

TITRES = [(0, "Apprenti·e calculateur·rice"), (40, "Arpenteur·euse de symétries"),
          (70, "Disciple d'Abel"), (90, "Héritier·ère de Galois")]


def _lire(entree):
    try:
        return entree("> ")
    except EOFError:
        return "q"


def generer(nb, graine=None):
    rng = random.Random(graine)
    familles = QUESTIONS * (nb // len(QUESTIONS) + 1)
    rng.shuffle(familles)
    return [f(rng) for f in familles[:nb]]


def jouer(nb=8, graine=None, entree=input, sortie=print):
    """Lance une partie.  ``entree``/``sortie`` sont injectables pour les tests."""
    sortie("\n🎲  LE DÉFI DE GALOIS  🎲")
    sortie("Tape ta réponse, « ? » pour un indice (−½ point) ou « q » pour abandonner.\n")
    score, serie, meilleure_serie = 0.0, 0, 0
    questions = generer(nb, graine)
    for num, q in enumerate(questions, 1):
        sortie(f"Question {num}/{nb} — {q.enonce}")
        points = 1.0
        rep = _lire(entree)
        if rep.strip() == "?":
            sortie(f"   💡 Indice : {q.indice}")
            points = 0.5
            rep = _lire(entree)
        if rep.strip().lower() == "q":
            sortie("Abandon... Galois, lui, a rédigé son testament mathématique en une seule nuit !")
            break
        if _normaliser(rep) == _normaliser(q.reponse):
            score += points
            serie += 1
            meilleure_serie = max(meilleure_serie, serie)
            bonus = "  🔥 série de " + str(serie) if serie >= 3 else ""
            sortie(f"   ✔ Bravo !{bonus}  {q.explication}\n")
        else:
            serie = 0
            sortie(f"   ✘ Raté, la réponse était « {q.reponse} ». {q.explication}\n")
    pourcentage = round(100 * score / nb)
    titre = [t for seuil, t in TITRES if pourcentage >= seuil][-1]
    sortie(f"Score final : {score:g}/{nb} ({pourcentage} %) — meilleure série : {meilleure_serie}")
    sortie(f"Ton titre : 🏅 {titre}")
    return score
