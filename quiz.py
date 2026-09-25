"""Le Défi de Galois : un quiz à questions tirées au hasard.

Chaque question est calculée par la bibliothèque elle-même : la réponse
attendue n'est jamais écrite à la main, elle est *démontrée* par le programme.
"""

import random
from dataclasses import dataclass

from catalogue import (cyclique, depuis_cycles, diedral, horloge, inversibles_mod,
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


QUESTIONS = [q_horloge, q_inverse_modulo, q_ordre_permutation, q_commutatif,
             q_cardinal_sn, q_symetries_polygone, q_composition_polygone,
             q_melange, q_lagrange]

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
