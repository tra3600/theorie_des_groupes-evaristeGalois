"""La théorie des groupes, d'Évariste Galois à ton écran.

Un petit laboratoire interactif : des cas illustrés (horloge, triangle,
carré, cartes à jouer, faux groupes, quaternions, équation du 5e degré)
et un quiz pour s'entraîner.

Exemples :
    python groups-theory-EGalois.py                    # menu interactif
    python groups-theory-EGalois.py --cas triangle     # un cas précis
    python groups-theory-EGalois.py --tout --sauver figures
    python groups-theory-EGalois.py --quiz 10
    python groups-theory-EGalois.py --tout --sans-graphique
"""

import argparse
import os
import sys

import catalogue as cat
from catalogue import depuis_cycles, notation_cyclique

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")


# ----------------------------------------------------------------------
# Affichage
# ----------------------------------------------------------------------
class Afficheur:
    """Gère l'affichage des figures : écran, fichier PNG, ou rien du tout."""

    def __init__(self, graphique=True, dossier=None):
        self.graphique = graphique or dossier is not None
        self.dossier = dossier
        if dossier:
            import matplotlib
            matplotlib.use("Agg")
            os.makedirs(dossier, exist_ok=True)

    def __call__(self, fabrique, nom):
        if not self.graphique:
            return
        import matplotlib.pyplot as plt
        fig = fabrique()
        if self.dossier:
            chemin = os.path.join(self.dossier, f"{nom}.png")
            fig.savefig(chemin, dpi=130)
            print(f"   🖼  figure enregistrée : {chemin}")
            plt.close(fig)
        else:
            plt.show()


def titre(texte):
    print("\n" + "═" * 72)
    print(f"  {texte}")
    print("═" * 72)


def para(texte):
    print(texte.strip("\n"))


# ----------------------------------------------------------------------
# Les cas illustrés
# ----------------------------------------------------------------------
def cas_horloge(montrer):
    import illustrations as ill
    titre("1. L'HORLOGE : ton premier groupe sans le savoir")
    G = cat.horloge()
    para("""
Il est 9h. Dans 5 heures, il sera... 14h ? Non : 2h ! Sur un cadran, on compte
« modulo 12 ». Les 12 heures munies de l'addition forment le groupe Z/12Z.
""")
    for a, b in [(9, 5), (11, 13), (3, -7)]:
        r = G.op(a % 12, b % 12)
        print(f"   {a}h {'+' if b >= 0 else '−'} {abs(b)}h  →  {r or 12}h")
    print("\nLes 4 règles du jeu sont-elles respectées ?")
    print(G.verifier())
    print(f"\nCombien d'heures pour faire un tour complet en avançant de 3h en 3h ? "
          f"{G.ordre_element(3)} sauts.")
    print(f"Et de 5h en 5h ? {G.ordre_element(5)} sauts : 5 visite TOUTES les heures, "
          f"c'est un générateur (car pgcd(5, 12) = 1).")
    gens = [a for a in G.elements if G.ordre_element(a) == 12]
    print(f"Les générateurs de l'horloge : {', '.join(G.etiquette(a) for a in gens)}")
    montrer(lambda: ill.horloge(9, 17), "01_horloge")
    montrer(lambda: ill.graphe_cayley(G, [1, 5]), "01_horloge_cayley")


def cas_triangle(montrer):
    import illustrations as ill
    titre("2. LE TRIANGLE : symétries et permutations (D3 ≅ S3)")
    D3, S3 = cat.diedral(3), cat.symetrique(3)
    para("""
Découpe un triangle équilatéral dans du carton : de combien de façons peux-tu le
reposer dans son trou ? 3 rotations (0°, 120°, 240°) et 3 retournements (miroirs).
Chaque symétrie mélange les sommets 1, 2, 3 : c'est une permutation !
""")
    for p in D3:
        print(f"   {D3.etiquette(p):>6}  mélange les sommets comme  {notation_cyclique(p)}")
    print("\nTable de Cayley de D3 :")
    print(D3.table_texte())
    ce = D3.contre_exemple_commutativite()
    a, b = ce
    print(f"\n⚠  L'ordre compte : {D3.etiquette(a)}∘{D3.etiquette(b)} = {D3.etiquette(D3.op(a, b))}"
          f" mais {D3.etiquette(b)}∘{D3.etiquette(a)} = {D3.etiquette(D3.op(b, a))}.")
    print("   Retourner puis tourner ≠ tourner puis retourner : le groupe n'est PAS commutatif.")
    iso = D3.isomorphisme(S3)
    print(f"\nD3 et S3 ont-ils la même structure ? {'OUI' if iso else 'non'} — "
          "les 6 symétries du triangle SONT les 6 permutations de 3 objets.")
    print("C'est exactement le groupe de Galois de l'équation x³ = 2 (voir le cas « galois »).")
    montrer(lambda: ill.symetries_polygone(3), "02_triangle_symetries")
    montrer(lambda: ill.composition_polygone(3, D3.elements[1], D3.elements[3]), "02_triangle_composition")
    montrer(lambda: ill.table_cayley(D3), "02_triangle_table")


def cas_carre(montrer):
    import illustrations as ill
    titre("3. LE CARRÉ : sous-groupes, Lagrange et sous-groupes distingués")
    G = cat.diedral(4)
    para("""
Le carré a 8 symétries : 4 rotations et 4 miroirs (2 diagonales, 2 médianes).
Certaines parties de D4 forment un groupe à elles seules : ce sont des sous-groupes.
""")
    sgs = G.sous_groupes()
    for H in sgs:
        noms = "{" + ", ".join(G.etiquette(x) for x in H) + "}"
        marque = "  ← distingué" if G.est_distingue(H) else ""
        print(f"   ordre {len(H)} : {noms}{marque}")
    tailles = sorted({len(H) for H in sgs})
    print(f"\n🔎 Théorème de Lagrange : les tailles possibles sont {tailles}, "
          f"toutes des diviseurs de {len(G)}.")
    rotations = G.engendre(G.elements[1])
    print("\nPourquoi ? Les classes aH découpent le groupe en paquets de même taille :")
    for c in G.classes_a_gauche(rotations):
        print("   {" + ", ".join(G.etiquette(x) for x in c) + "}")
    print("\nClasses de conjugaison (les éléments « qui se ressemblent ») :")
    for c in G.classes_de_conjugaison():
        print("   {" + ", ".join(G.etiquette(x) for x in c) + "}")
    montrer(lambda: ill.symetries_polygone(4), "03_carre_symetries")
    montrer(lambda: ill.treillis_sous_groupes(G), "03_carre_sous_groupes")
    montrer(lambda: ill.graphe_cayley(G, [G.elements[1], G.elements[4]], "anneaux"),
            "03_carre_cayley")


def cas_cartes(montrer):
    import illustrations as ill
    titre("4. LES CARTES : le tour de magie du mélange parfait")
    para("""
Un « faro » parfait : on coupe le paquet en deux moitiés égales et on intercale
les cartes une à une. Les magiciens savent qu'avec 52 cartes, 8 mélanges parfaits
ramènent le paquet exactement dans son ordre de départ. Magie ? Non : théorie des groupes !
""")
    for n in (8, 10, 52):
        for interieur, nom in ((False, "extérieur"), (True, "intérieur")):
            p = cat.melange_parfait(n, interieur)
            k = len(cat.groupe_engendre_par_permutation(p))
            print(f"   {n:>2} cartes, faro {nom:<9} : retour après {k:>2} mélanges")
    p = cat.melange_parfait(8)
    print(f"\nPour 8 cartes, le mélange extérieur est la permutation {notation_cyclique(p)} :")
    longueurs = sorted(len(c) for c in cat.cycles(p))
    print(f"   ses cycles ont pour longueurs {longueurs} (+ {n_fixes(p)} cartes immobiles)"
          f" → ordre = ppcm{tuple(longueurs)} = {len(cat.groupe_engendre_par_permutation(p))}.")
    montrer(lambda: ill.melanges(cat.melange_parfait(52), "Faro extérieur"), "04_cartes_52")
    montrer(lambda: ill.melanges(cat.melange_parfait(52, True), "Faro intérieur"), "04_cartes_52_interieur")


def cas_contre_exemples(montrer):
    titre("5. LES IMPOSTEURS : quand ce n'est PAS un groupe")
    para("""
Un groupe doit respecter 4 règles : fermeture, associativité, élément neutre,
inverses. Voici des candidats qui trichent — le programme trouve la faille !
""")
    for G in (cat.faux_addition_bornee(), cat.faux_soustraction(4), cat.faux_maximum(),
              cat.faux_multiplication_avec_zero(5), cat.faux_multiplication_mod(6),
              cat.inversibles_mod(6), cat.faux_multiplication_mod(7)):
        print(f"\n▶ {G.nom}")
        print(G.verifier())
    para("""
Morale : (Z/nZ \\ {0}, ×) est un groupe exactement quand n est premier.
Pour n = 6, il faut ne garder que les nombres premiers avec 6 : (Z/6Z)* = {1, 5}.
""")


def cas_quatre(montrer):
    import illustrations as ill
    titre("6. MÊME TAILLE, AUTRE GROUPE : Z/4Z, Klein et les quaternions")
    Z4, V4, U8 = cat.cyclique(4), cat.klein(), cat.inversibles_mod(8)
    para("""
Deux groupes à 4 éléments sont-ils forcément « pareils » ? Comptons les éléments
de chaque ordre : c'est leur empreinte digitale.
""")
    for G in (Z4, V4, U8):
        print(f"   {G.nom:<8} ordres → {G.signature()}   cyclique : {'oui' if G.est_cyclique() else 'non'}")
    print("\n→ Z/4Z a un élément d'ordre 4, pas V4 : ils sont DIFFÉRENTS.")
    print(f"→ (Z/8Z)* ≅ V4 ? {'oui' if U8.isomorphisme(V4) else 'non'} : 3² ≡ 5² ≡ 7² ≡ 1 (mod 8).")
    Q8, D4 = cat.quaternions(), cat.diedral(4)
    para("""
À 8 éléments, deux groupes non commutatifs : D4 (le carré) et Q8, les quaternions
que Hamilton grava sur un pont de Dublin en 1843 : i² = j² = k² = ijk = −1.
""")
    print(Q8.table_texte())
    print(f"\n   D4 ordres → {D4.signature()}")
    print(f"   Q8 ordres → {Q8.signature()}")
    print("→ D4 a 5 éléments d'ordre 2, Q8 un seul (−1) : ce ne sont pas les mêmes groupes.")
    print(f"→ Dans Q8, TOUS les {len(Q8.sous_groupes())} sous-groupes sont distingués, "
          "bien que Q8 ne soit pas commutatif. Curieux, non ?")
    montrer(lambda: ill.table_cayley(Q8), "06_quaternions_table")
    montrer(lambda: ill.graphe_cayley(Q8, [(0, 1, 0, 0), (0, 0, 1, 0)], "anneaux"),
            "06_quaternions_cayley")


def cas_galois(montrer):
    import illustrations as ill
    titre("7. GALOIS ET L'ÉQUATION DU 5e DEGRÉ")
    para("""
Pour x² + bx + c = 0, on connaît la formule avec √. Pour les degrés 3 et 4 aussi
(Cardan, Ferrari, XVIe siècle). Et le degré 5 ? Pendant 300 ans, personne n'a trouvé.

À 20 ans, dans la nuit précédant son duel fatal (29-30 mai 1832), Évariste Galois explique
pourquoi : à chaque équation on associe un groupe de permutations de ses racines.
L'équation se résout avec des racines carrées, cubiques... si et seulement si
ce groupe est « résoluble » : en prenant les commutateurs aba⁻¹b⁻¹ encore et
encore, on doit finir par tomber sur {e}.

Exemple : les 3 racines de x³ = 2 sont ∛2, ∛2·ω, ∛2·ω² (ω = e^(2iπ/3)).
Leur groupe de Galois est S3 tout entier.
""")
    groupes = [cat.symetrique(3), cat.symetrique(4), cat.symetrique(5)]
    for G in groupes:
        serie = G.serie_derivee()
        chaine = " ⊃ ".join(str(len(H)) for H in serie)
        verdict = "résoluble ✔" if len(serie[-1]) == 1 else "NON résoluble ✘"
        print(f"   {G.nom} : tailles {chaine}  →  {verdict}")
    A5 = cat.alterne(5)
    tailles = sorted(len(c) for c in A5.classes_de_conjugaison())
    para(f"""
Pourquoi S5 coince-t-il ? Son dérivé est A5 (60 éléments), et le dérivé de A5 est
A5 lui-même. Mieux : A5 est « simple ». Preuve express par le programme : ses
classes de conjugaison ont pour tailles {tailles}. Un sous-groupe distingué est une
réunion de classes contenant {{e}} dont la taille divise 60 — et seules les
réunions de taille 1 et 60 conviennent. Pas d'échelle pour descendre !
""")
    combinaisons = _sommes_de_classes(tailles)
    print(f"   Tailles possibles (1 + réunions de classes) : {sorted(combinaisons)}")
    print(f"   Celles qui divisent 60 : {sorted(t for t in combinaisons if 60 % t == 0)}")
    x = depuis_cycles(5, (1, 2, 3, 4, 5))
    print(f"\nExemple concret : x⁵ − x − 1 = 0 a pour groupe de Galois S5 tout entier.")
    print(f"   (il contient le 5-cycle {notation_cyclique(x)} et une transposition, qui engendrent S5 :"
          f" {len(cat.symetrique(5).engendre(x, depuis_cycles(5, (1, 2))))} éléments)")
    print("   ⇒ AUCUNE formule avec des radicaux ne donne ses racines. Merci Galois !")
    montrer(lambda: ill.series_derivees(groupes), "07_galois_series")
    montrer(lambda: ill.treillis_sous_groupes(cat.symetrique(3)), "07_galois_s3_sous_groupes")


def n_fixes(p):
    return sum(1 for i, x in enumerate(p) if i == x)


def _sommes_de_classes(tailles):
    autres = list(tailles)
    autres.remove(1)
    sommes = {1}
    for t in autres:
        sommes |= {s + t for s in sommes}
    return sommes


CAS = {
    "horloge": ("L'horloge : Z/12Z", cas_horloge),
    "triangle": ("Les symétries du triangle : D3 ≅ S3", cas_triangle),
    "carre": ("Le carré : sous-groupes et Lagrange", cas_carre),
    "cartes": ("Le tour de magie du mélange parfait", cas_cartes),
    "imposteurs": ("Les faux groupes et leurs failles", cas_contre_exemples),
    "quatre": ("Même taille, autre groupe : V4, Q8...", cas_quatre),
    "galois": ("Galois et l'équation du 5e degré", cas_galois),
}


def menu(montrer):
    import quiz
    while True:
        print("\n╔══════════════════════════════════════════════╗")
        print("║   🎓  LE LABORATOIRE D'ÉVARISTE GALOIS  🎓   ║")
        print("╚══════════════════════════════════════════════╝")
        for i, (cle, (desc, _)) in enumerate(CAS.items(), 1):
            print(f"  {i}. {desc}")
        print("  Q. 🎲 Le Défi de Galois (quiz)")
        print("  T. Tout voir d'un coup")
        print("  X. Quitter")
        try:
            choix = input("Ton choix : ").strip().lower()
        except EOFError:
            choix = "x"
        cles = list(CAS)
        if choix in ("x", "quitter", ""):
            print("À bientôt ! « Je n'ai pas le temps » — É. Galois, 29 mai 1832.")
            return
        if choix == "q":
            quiz.jouer()
        elif choix == "t":
            for _, f in CAS.values():
                f(montrer)
        elif choix.isdigit() and 1 <= int(choix) <= len(cles):
            CAS[cles[int(choix) - 1]][1](montrer)
        elif choix in CAS:
            CAS[choix][1](montrer)
        else:
            print("Choix inconnu.")


def main(argv=None):
    parser = argparse.ArgumentParser(description="Théorie des groupes ludique et illustrée.")
    parser.add_argument("--cas", choices=list(CAS), action="append",
                        help="lancer un cas illustré (répétable)")
    parser.add_argument("--tout", action="store_true", help="lancer tous les cas")
    parser.add_argument("--quiz", type=int, nargs="?", const=8, metavar="N",
                        help="jouer au quiz (N questions, 8 par défaut)")
    parser.add_argument("--graine", type=int, help="graine du hasard pour le quiz")
    parser.add_argument("--sans-graphique", action="store_true", help="texte seulement")
    parser.add_argument("--sauver", metavar="DOSSIER", help="enregistrer les figures en PNG")
    args = parser.parse_args(argv)

    montrer = Afficheur(graphique=not args.sans_graphique, dossier=args.sauver)
    cas = list(CAS) if args.tout else (args.cas or [])
    for c in cas:
        CAS[c][1](montrer)
    if args.quiz:
        import quiz
        quiz.jouer(args.quiz, graine=args.graine)
    if not cas and not args.quiz:
        menu(montrer)


if __name__ == "__main__":
    main()
