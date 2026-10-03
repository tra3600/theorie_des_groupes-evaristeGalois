"""Les cas avancés du laboratoire : quotients, classification, Sylow, Cayley, Burnside,
cryptographie, Rubik et groupes de Galois de polynômes.

Chaque cas reçoit ``montrer(fabrique, nom)`` et affiche des faits *calculés* par le programme.
"""

import math

import arithmetique as ar
import catalogue as cat
import denombrement as dn
import galois_polynomes as gp
import permgroup as pg
from catalogue import notation_cyclique


def titre(texte):
    print("\n" + "═" * 72)
    print(f"  {texte}")
    print("═" * 72)


def para(texte):
    print(texte.strip("\n"))


def liste(G, H):
    return "{" + ", ".join(G.etiquette(x) for x in H) + "}"


def espace(n):
    """12345678 -> '12 345 678'."""
    return f"{n:,}".replace(",", " ")


def _facteurs(n):
    return " · ".join(f"{p}^{e}" if e > 1 else str(p) for p, e in sorted(ar.facteurs_premiers(n).items()))


# ----------------------------------------------------------------------
def cas_quotients(montrer):
    import illustrations as ill
    titre("8. CENTRE, QUOTIENTS ET HOMOMORPHISMES")
    para("""
Le centre Z(G) rassemble les éléments qui commutent avec tout le monde. Quand un
sous-groupe N est distingué, on peut « recoller » ses classes en un nouveau groupe,
le quotient G/N. Un homomorphisme est une fonction qui respecte l'opération.
""")
    D4, Q8 = cat.diedral(4), cat.quaternions()
    for G in (D4, Q8, cat.symetrique(3)):
        z, classes = G.equation_des_classes()
        somme = " + ".join(map(str, classes))
        print(f"   {G.nom:<3} centre {liste(G, G.centre())}\n"
              f"       équation des classes : {len(G)} = {z}" + (f" + {somme}" if classes else "")
              + "   (chaque classe a un cardinal qui divise |G|)")
    print("\nQuotients (le résultat est un nouveau groupe, comparé à un groupe connu) :")
    V4 = cat.klein()
    S4 = cat.symetrique(4)
    klein_s4 = [x for x in S4 if x in (cat.depuis_cycles(4), cat.depuis_cycles(4, (1, 2), (3, 4)),
                                       cat.depuis_cycles(4, (1, 3), (2, 4)), cat.depuis_cycles(4, (1, 4), (2, 3)))]
    A4 = cat.alterne(4)
    cas = [(D4, D4.centre(), "D4/Z(D4)", V4, "V4"), (Q8, Q8.centre(), "Q8/{±1}", V4, "V4"),
           (S4, klein_s4, "S4/V4", cat.symetrique(3), "S3"),
           (A4, klein_s4, "A4/V4", cat.cyclique(3), "Z/3Z"), (cat.cyclique(12), cat.cyclique(12).engendre(4),
                                                              "Z/12Z/<4>", cat.cyclique(4), "Z/4Z")]
    for G, N, nom, but, nom_but in cas:
        Q = G.quotient(N, nom)
        print(f"   {nom:<10} a {len(Q)} éléments ({len(G)}/{len(N)}) et ≅ {nom_but} : "
              f"{'OUI' if Q.isomorphisme(but) else 'non'}")
    S3, Z2 = cat.symetrique(3), cat.cyclique(2)

    def signe(p):
        return 0 if cat.signature_permutation(p) == 1 else 1
    ker = S4.noyau(signe, Z2)
    print(f"\nHomomorphisme « signature » S4 → Z/2Z : est-ce un homomorphisme ? {S4.est_homomorphisme(signe, Z2)}")
    print(f"   noyau = A4 ({len(ker)} éléments), image de taille {len(S4.image(signe))} : "
          f"{len(S4)} = {len(ker)} × {len(S4.image(signe))}  (premier théorème d'isomorphisme)")
    Z12 = cat.cyclique(12)
    print(f"   x ↦ x mod 4 de Z/12Z vers Z/4Z : homomorphisme ? {Z12.est_homomorphisme(lambda x: x % 4, cat.cyclique(4))}, "
          f"noyau {liste(Z12, Z12.noyau(lambda x: x % 4, cat.cyclique(4)))}")
    print(f"   x ↦ x² de Z/12Z vers Z/12Z : homomorphisme ? {Z12.est_homomorphisme(lambda x: x * x % 12, Z12)}")
    print("\nGroupes simples (aucun quotient non trivial : les « atomes » de la théorie) :")
    for G in (cat.cyclique(7), cat.cyclique(6), cat.diedral(5), A4, S4, cat.alterne(5)):
        print(f"   {G.nom:<5} simple ? {'OUI' if G.est_simple() else 'non'}")
    print("   Les groupes abéliens simples sont les Z/pZ ; le premier simple non abélien est A5 (60 éléments).")
    montrer(lambda: ill.table_cayley(D4.quotient(D4.centre(), "D4/Z")), "08_quotient_d4")


def cas_classification(montrer):
    import illustrations as ill
    titre("9. TOUS LES GROUPES D'ORDRE ≤ 12, ET LEUR PETITE CARTE D'IDENTITÉ")
    para("""
Combien existe-t-il de groupes à n éléments, à isomorphisme près ? Le programme les
fabrique, vérifie qu'ils sont bien des groupes, et qu'aucun n'est isomorphe à un autre.
""")
    attendu = [1, 1, 1, 2, 1, 2, 1, 5, 2, 2, 1, 5]
    par_ordre = {n: cat.groupes_d_ordre(n) for n in range(1, 13)}
    for n, groupes in par_ordre.items():
        assert all(len(G) == n and G.verifier().est_groupe for G in groupes)
        for i, G in enumerate(groupes):
            assert all(G.isomorphisme(H) is None for H in groupes[i + 1:])
    print(f"   {'ordre':>5} {'nombre':>7}  groupes")
    for n, groupes in par_ordre.items():
        print(f"   {n:>5} {len(groupes):>7}  {', '.join(G.nom for G in groupes)}")
    total = sum(len(g) for g in par_ordre.values())
    print(f"\n   Total : {total} groupes, comme attendu ({attendu}) : "
          f"{'✔' if [len(g) for g in par_ordre.values()] == attendu else '✘'}")
    print("\nTrois leçons :")
    print("   • Pour un ordre premier, un seul groupe : Z/pZ (Lagrange : tout élément non neutre engendre).")
    Z2, Z3, Z4 = cat.cyclique(2), cat.cyclique(3), cat.cyclique(4)
    print(f"   • Z/6Z ≅ Z/2Z × Z/3Z ? {'OUI' if cat.cyclique(6).isomorphisme(Z2.produit_direct(Z3)) else 'non'}"
          f" (théorème chinois : 2 et 3 sont premiers entre eux)")
    print(f"     Z/4Z ≅ Z/2Z × Z/2Z ? {'OUI' if Z4.isomorphisme(Z2.produit_direct(Z2)) else 'non'}"
          " (2 et 2 ne sont pas premiers entre eux)")
    print(f"   • À l'ordre 6, la non-commutativité apparaît : S3 est le plus petit groupe non abélien.")
    print("\nEt l'ordre 24 ? Il y a 15 groupes. Même taille, empreintes différentes :")
    for G in (cat.symetrique(4), cat.sl2_mod(3), cat.alterne(4).produit_direct(Z2), cat.cyclique(24)):
        print(f"   {G.nom:<11} ordres des éléments {G.signature()}  résoluble : {'oui' if G.est_resoluble() else 'non'}")
    print("   (S4 et SL(2,3) ont chacun 24 éléments mais SL(2,3) n'a qu'une seule symétrie d'ordre 2.)")
    montrer(lambda: ill.classification(par_ordre), "09_classification")


def cas_sylow(montrer):
    import illustrations as ill
    titre("10. SYLOW : COMBIEN DE SOUS-GROUPES DE TAILLE p^a ?")
    para("""
Lagrange dit qu'une taille de sous-groupe divise |G|, mais il ne dit pas que chaque
diviseur est atteint. Sylow (1872) comble une partie du vide : si p^a est la plus grande
puissance de p qui divise |G|, il EXISTE un sous-groupe de taille p^a, et leur nombre n_p
vérifie  n_p ≡ 1 (mod p)  et  n_p divise |G|/p^a. Le programme les construit.
""")
    print(f"   {'groupe':<6}{'|G|':>5}{'p':>4}{'p^a':>5}{'n_p':>6}  n_p ≡ 1 (mod p) ?  n_p | |G|/p^a ?")
    for G in (cat.symetrique(3), cat.diedral(4), cat.alterne(4), cat.symetrique(4), cat.diedral(6),
              cat.alterne(5), cat.symetrique(5)):
        n = len(G)
        for p in ar.facteurs_premiers(n):
            P, nb = G.sylow(p)
            reste = n // len(P)
            assert nb % p == 1 and reste % nb == 0
            print(f"   {G.nom:<6}{n:>5}{p:>4}{len(P):>5}{nb:>6}  {'✔':>13}  {'✔':>18}")
    A4 = cat.alterne(4)
    tailles = sorted({len(H) for H in A4.sous_groupes()})
    print(f"\n   Le sous-groupe de Sylow 2 de A4 est unique (n_2 = 1) : donc distingué. C'est V4.")
    print(f"   Tailles des sous-groupes de A4 : {tailles}. Il divise 12 mais n'a AUCUN sous-groupe d'ordre 6 :")
    print("   la réciproque de Lagrange est fausse !")
    print("\n   Application : un groupe d'ordre 15 est-il forcément cyclique ?")
    print("      n_5 ≡ 1 (mod 5) et n_5 | 3  ⇒ n_5 = 1 ;  n_3 ≡ 1 (mod 3) et n_3 | 5  ⇒ n_3 = 1.")
    print("      Les deux Sylow sont distingués, leur produit est direct : Z/3 × Z/5 ≅ Z/15. Un seul groupe d'ordre 15.")
    montrer(lambda: ill.treillis_sous_groupes(A4), "10_sylow_a4")


def cas_cayley(montrer):
    titre("11. LE THÉORÈME DE CAYLEY : TOUT GROUPE EST UN GROUPE DE PERMUTATIONS")
    para("""
Cayley (1854) : chaque groupe fini G s'envoie dans S_|G| en associant à g la permutation
« multiplier à gauche par g ». C'est un homomorphisme injectif : G se retrouve, fidèlement,
parmi les permutations. Étudier les groupes, c'est donc étudier les permutations.
""")
    for G in (cat.cyclique(5), cat.symetrique(3), cat.diedral(4), cat.quaternions()):
        rep = G.representation_reguliere()
        elements = list(G.elements)
        ok = all(rep[G.op(a, b)] == cat.composer(rep[a], rep[b]) for a in elements for b in elements)
        injectif = len(set(rep.values())) == len(G)
        image = pg.GroupePermutations(list(rep.values()))
        print(f"   {G.nom:<4}: |G| = {len(G)}, dans S{len(G)} : homomorphisme {'✔' if ok else '✘'}, "
              f"injectif {'✔' if injectif else '✘'}, la copie a {image.ordre()} éléments")
    D4 = cat.diedral(4)
    r = D4.representation_reguliere()
    g = D4.elements[1]
    print(f"\n   Exemple : la rotation {D4.etiquette(g)} de D4 devient la permutation {notation_cyclique(r[g])} de 8 objets.")
    print(f"   Mais D4 tient déjà dans S4 (son action sur les 4 sommets) : {notation_cyclique(g)}.")
    print("   Cayley donne une plongée toujours valable, pas toujours la plus économe.")


def cas_burnside(montrer):
    import illustrations as ill
    titre("12. BURNSIDE : COMPTER LES COLORIAGES À SYMÉTRIE PRÈS")
    para("""
Combien de colliers de 6 perles en 2 couleurs ? 2⁶ = 64 si les perles sont numérotées, mais
une rotation redonne le même collier. Burnside : le nombre de coloriages distincts est la
moyenne, sur toutes les symétries g, du nombre de coloriages que g laisse inchangés,
soit k^(nombre de cycles de g).
""")
    print(f"   {'n':>3} {'k':>3} {'colliers (C_n)':>16} {'bracelets (D_n)':>17}  vérification par énumération")
    for n, k in ((4, 2), (5, 2), (6, 2), (6, 3), (7, 2)):
        rot = [cat.rotation(n, j) for j in range(n)]
        dih = list(cat.diedral(n).elements)
        c, b = dn.colliers(n, k), dn.bracelets(n, k)
        ok = dn.orbites_force_brute(rot, n, k) == c and dn.orbites_force_brute(dih, n, k) == b
        print(f"   {n:>3} {k:>3} {c:>16} {b:>17}  {'✔ identique' if ok else '✘'}")
    G = dn.groupe_des_rotations_du_cube()
    S4 = cat.symetrique(4)
    print(f"\nLes rotations du cube : {len(G)} éléments, ≅ S4 ? {'OUI' if G.isomorphisme(S4) else 'non'}"
          " (elles permutent les 4 grandes diagonales du cube).")
    perms = [pg.action_sur("faces", m) for m in pg.rotations_cube()[0]]
    print("   Types de cycles des rotations sur les 6 faces : " +
          ", ".join(f"{'+'.join(map(str, t))}×{c}" for t, c in dn.indice_des_cycles(perms).items()))
    print(f"\n   {'objet':<12}" + "".join(f"{k:>7} coul." for k in (2, 3, 4)))
    donnees = {}
    for genre, nom in (("faces", "faces"), ("sommets", "sommets"), ("aretes", "arêtes")):
        valeurs = {k: dn.coloriages_du_cube(genre, k) for k in (2, 3, 4)}
        donnees[nom] = valeurs
        print(f"   {nom:<12}" + "".join(f"{valeurs[k]:>12}" for k in (2, 3, 4)))
    force = dn.orbites_force_brute(perms, 6, 3)
    print(f"\n   Vérification : 3^6 = 729 coloriages de faces rangés en orbites par force brute → {force} "
          f"(Burnside : {dn.coloriages_du_cube('faces', 3)}).")
    montrer(lambda: ill.barres_burnside(donnees, "Colorier un cube à rotation près",
                                        "peindre faces, sommets ou arêtes : le lemme de Burnside compte les cubes différents"),
            "12_burnside_cube")


def cas_crypto(montrer):
    import illustrations as ill
    titre("13. LES GROUPES (Z/nZ)* DANS TON NAVIGATEUR : DIFFIE–HELLMAN ET RSA")
    para("""
(Z/nZ)* est le groupe des nombres premiers avec n, pour la multiplication. Lagrange donne
a^φ(n) ≡ 1 (mod n) (Euler) ; c'est le moteur de la cryptographie à clé publique.
""")
    p, g = 23, 5
    G = cat.inversibles_mod(p)
    print(f"   (Z/{p}Z)* a {len(G)} éléments ; cyclique ? {'oui' if G.est_cyclique() else 'non'} ; "
          f"{len(ar.racines_primitives(p))} générateurs = φ({p - 1}) = {ar.phi(p - 1)}.")
    a, b = 6, 15
    A, B, ka, kb = ar.echange_diffie_hellman(p, g, a, b)
    print(f"\n   Diffie–Hellman, p = {p}, g = {g} : Alice (secret {a}) envoie {A}, Bob (secret {b}) envoie {B}.")
    print(f"   Clé d'Alice = {B}^{a} = {ka} ; clé de Bob = {A}^{b} = {kb} → même clé, sans l'avoir transmise.")
    print("   Ève voit p, g, A, B et doit résoudre un logarithme discret : g^x ≡ A (mod p).")
    P = 1_000_003
    g2 = ar.plus_petite_racine_primitive(P)
    x = 654_321
    h = pow(g2, x, P)
    naif, essais = ar.logarithme_naif(g2, h, P)
    bsgs, ops = ar.logarithme_discret(g2, h, P)
    print(f"\n   Pour p = {espace(P)} (g = {g2}, secret {espace(x)}), Ève cherche x :")
    print(f"      en essayant toutes les puissances : x = {naif} après {espace(essais)} essais")
    print(f"      par pas de bébé / pas de géant     : x = {bsgs} après {espace(ops)} opérations")
    print(f"   Environ la racine carrée du nombre de cas ({espace(ops)} ≈ √{espace(P)}). "
          "Avec un p de 600 chiffres, même une racine carrée est hors d'atteinte.")
    print("\n   RSA avec p = 61, q = 53, e = 17, message 65 :")
    n, f, d, c, m = ar.rsa_jouet(61, 53, 17, 65)
    print(f"      n = {n}, φ(n) = {f}, d = e⁻¹ mod φ(n) = {d} ; chiffré 65^17 = {c} ; déchiffré {c}^{d} = {m}.")
    print(f"      Euler : 65^φ(n) mod n = {pow(65, f, n)} ; l'exposant du groupe (Carmichael) λ(n) = {ar.carmichael(n)}"
          f" divise φ(n) et donne aussi une clé : d' = {pow(17, -1, ar.carmichael(n))}.")
    cycl = [n_ for n_ in range(2, 31) if cat.inversibles_mod(n_).est_cyclique()]
    assert cycl == [n_ for n_ in range(2, 31) if ar.a_une_racine_primitive(n_)]
    print(f"\n   (Z/nZ)* est cyclique pour n = {cycl} (n ≤ 30) : exactement 1, 2, 4, p^k, 2p^k (p premier impair).")
    montrer(lambda: ill.puissances_modulaires(13, [2, 3]), "13_multiplication_mod_13")


def cas_rubik(montrer):
    import illustrations as ill
    titre("14. LE RUBIK'S CUBE : UN GROUPE DE 43 TRILLIONS D'ÉLÉMENTS")
    para("""
Les 6 quarts de tour de face engendrent un groupe de permutations des 48 autocollants mobiles.
On ne peut pas le lister : on lit son ordre sur la chaîne de stabilisateurs de Schreier–Sims
(un produit de tailles d'orbites), sans énumérer un seul élément.
""")
    R, U = pg.mouvement_cube("R"), pg.mouvement_cube("U")
    cube = pg.cube_rubik(3)
    ordre = cube.ordre()
    print(f"   Rubik 3×3×3 : {espace(ordre)} positions atteignables")
    print(f"      tailles d'orbites de la chaîne : {cube.tailles_orbites()}")
    formule = math.factorial(8) * 3 ** 8 * math.factorial(12) * 2 ** 12 // 12
    print(f"      = 8!·3⁸·12!·2¹² / 12 = {espace(formule)} {'✔' if formule == ordre else '✘'}"
          "  (orientations des coins mod 3, des arêtes mod 2, parités des permutations liées)")
    print(f"      décomposition : {_facteurs(ordre)}")
    secondes = ordre
    annees = secondes / (365.25 * 24 * 3600)
    print(f"      une position par seconde : {annees:.2e} ans (l'univers a 1,4×10¹⁰ ans)")
    petit = pg.cube_rubik(2)
    print(f"\n   Rubik 2×2×2 (R, U, F) : {espace(petit.ordre())} positions = 7!·3⁶ = {espace(math.factorial(7) * 3 ** 6)}")
    print(f"   Rubik 3×3×3 avec seulement R et U : {espace(pg.GroupePermutations([R, U]).ordre())} positions")
    print("\n   Ordres de quelques algorithmes (nombre de répétitions pour revenir au départ) :")
    Ri, Ui = pg.inverse(R), pg.inverse(U)
    for nom, p in (("R", R), ("R U", pg.compose(R, U)), ("R U R' U'", pg.compose(pg.compose(R, U), pg.compose(Ri, Ui))),
                   ("R²", pg.compose(R, R))):
        print(f"      {nom:<10} ordre {pg.ordre_permutation(p)}")
    s = (pg.compose(R, U))
    x = pg.identite(len(R))
    for _ in range(105):
        x = pg.compose(s, x)
    print(f"      (R U)^105 = identité ? {x == pg.identite(len(R))}  : répéter « R U » 105 fois ramène le cube.")
    dist = cube.distribution_des_ordres(400)
    print(f"\n   Ordres de 400 positions tirées au hasard : de {min(dist)} à {max(dist)} (maximum théorique : 1260).")
    orbites = {tuple(cube.orbite(i)) for i in range(54)}
    print(f"   Les 54 autocollants se répartissent en {len(orbites)} orbites de tailles "
          f"{sorted(len(o) for o in orbites)} : 6 centres immobiles, puis les autocollants des coins "
          "et ceux des arêtes, qui ne se mélangent jamais entre eux.")
    montrer(lambda: ill.distribution_des_ordres(dist, "Ordre d'un mélange tiré au hasard sur le Rubik's Cube",
                                                "chaque barre : proportion de positions dont l'ordre est la valeur indiquée"),
            "14_rubik_ordres")


def cas_polynomes(montrer):
    titre("15. LE GROUPE DE GALOIS D'UN POLYNÔME, CALCULÉ (degrés 2 à 4)")
    para("""
Pour un polynôme à coefficients entiers, le discriminant dit si le groupe tient dans A_n ;
au degré 4, la « résolvante cubique » départage les cinq possibilités. Ces critères
(Kappe et Warren, 1989) sont programmés : on lit le groupe sur les coefficients.
""")
    exemples = [
        ("x² − 2", [1, 0, -2]), ("x² − 1", [1, 0, -1]),
        ("x³ − 2", [1, 0, 0, -2]), ("x³ − x − 1", [1, 0, -1, -1]), ("x³ − 3x + 1", [1, 0, -3, 1]),
        ("x³ − 1", [1, 0, 0, -1]),
        ("x⁴ − 2", [1, 0, 0, 0, -2]), ("x⁴ + 1", [1, 0, 0, 0, 1]),
        ("x⁴ + x³ + x² + x + 1", [1, 1, 1, 1, 1]), ("x⁴ + 8x + 12", [1, 0, 0, 8, 12]),
        ("x⁴ − x − 1", [1, 0, 0, -1, -1]), ("(x² − 2)(x² − 3)", [1, 0, -5, 0, 6]),
    ]
    print(f"   {'polynôme':<24}{'discriminant':>14}  {'groupe de Galois':<18}{'ordre':>6}")
    for nom, f in exemples:
        r = gp.groupe_de_galois(f)
        print(f"   {nom:<24}{gp.discriminant(f):>14}  {r['nom']:<18}{r['ordre']:>6}")
    print("\nPourquoi x³ − 2 a un groupe de 6 éléments : ses racines sont ∛2, ω∛2, ω²∛2,")
    r = gp.groupe_de_galois([1, 0, 0, -2])
    print("   et le programme répond : " + " ; ".join(r["raisons"]) + ".")
    print("\nQuelques raisonnements complets :")
    for nom, f in (("x⁴ − 2", [1, 0, 0, 0, -2]), ("x⁴ + x³ + x² + x + 1", [1, 1, 1, 1, 1]), ("x⁴ + 8x + 12", [1, 0, 0, 8, 12])):
        r = gp.groupe_de_galois(f)
        print(f"   {nom} → {r['nom']} :")
        for ligne in r["raisons"]:
            print(f"        · {ligne}")
    print("\nTous ces groupes sont des sous-groupes de S4, qui est résoluble : d'où les formules de Cardan et Ferrari.")
    for nom in ("S3", "D4", "A4", "S4"):
        G = gp.groupe_catalogue(nom)
        print(f"   {nom:<3} série dérivée : " + " ⊃ ".join(str(len(H)) for H in G.serie_derivee()))


def cas_quintique(montrer):
    import illustrations as ill
    titre("16. PROUVER QUE x⁵ − x − 1 N'EST PAS RÉSOLUBLE PAR RADICAUX")
    para("""
Au degré 5, plus de formule. Mais peut-on le *démontrer* pour un polynôme précis ? Dedekind :
si p ne divise pas le discriminant, la factorisation de f modulo p donne le type de cycles
d'une permutation du groupe de Galois. Il suffit de trouver
  • un premier où f reste irréductible (un 5-cycle : le groupe est transitif d'ordre divisible par 5),
  • un premier où f = (quadratique irréductible) × (3 linéaires) (une transposition),
car un sous-groupe transitif de S5 contenant une transposition est S5 tout entier.
""")
    f = [1, 0, 0, 0, -1, -1]
    D = gp.discriminant(f)
    print(f"   f = x⁵ − x − 1, discriminant {D} = {_facteurs(D)} (ces deux premiers sont exclus)")
    print(f"\n   {'p':>5}  degrés des facteurs de f mod p  → type de cycle d'un élément du groupe")
    for p in gp.premiers_jusqu_a(40):
        if D % p == 0:
            print(f"   {p:>5}  (p divise le discriminant)")
            continue
        t = gp.types_de_cycles_mod_p(f, p)
        print(f"   {p:>5}  {'+'.join(map(str, t)):<12}")
    ok, p5, ptr = gp.prouve_s5(f)
    print(f"\n   Preuve : irréductible modulo {p5} (un 5-cycle) ; modulo {ptr}, types {gp.types_de_cycles_mod_p(f, ptr)} "
          f"(une transposition) ⇒ Gal(f) = S5 : {ok}.")
    S5 = cat.symetrique(5)
    print(f"   S5 n'est pas résoluble (série dérivée {' ⊃ '.join(str(len(H)) for H in S5.serie_derivee())}) : "
          "aucune formule avec radicaux n'exprime les racines de f.")
    obs_a, tot_a = gp.statistique_de_frobenius(f, 3000)
    th_a = gp.frequences_theoriques(list(S5))
    g2 = [1, 0, 0, 0, 0, -2]
    obs_b, tot_b = gp.statistique_de_frobenius(g2, 3000)
    F20 = cat.permutations_du_groupe_affine(5)
    th_b = gp.frequences_theoriques(F20)
    print(f"\n   Chebotarev : avec {tot_a} premiers (p ≤ 3000), les types observés suivent les proportions du groupe :")
    print(f"   {'type':<12}{'observé S5':>12}{'théorie S5':>12}")
    for t in sorted(th_a, key=lambda t: (-len(t), t), reverse=True):
        print(f"   {'+'.join(map(str, t)):<12}{100 * obs_a.get(t, 0) / tot_a:>11.1f}%{100 * th_a[t]:>11.1f}%")
    G20 = cat.groupe_affine(5)
    print(f"\n   Contraste avec x⁵ − 2 : jamais de transposition parmi {tot_b} premiers. Son groupe est AGL(1,5), "
          f"d'ordre {len(G20)}, résoluble :")
    print(f"   série dérivée " + " ⊃ ".join(str(len(H)) for H in G20.serie_derivee()) +
          f" ; ses types de cycles : {sorted('+'.join(map(str, t)) for t in th_b)}.")
    print("   Ses racines s'écrivent avec des radicaux (⁵√2 et les racines 5-ièmes de l'unité) — celles de x⁵ − x − 1, jamais.")
    montrer(lambda: ill.frobenius({"x⁵ − x − 1 (groupe S5)": (obs_a, tot_a, th_a),
                                   "x⁵ − 2 (groupe AGL(1,5))": (obs_b, tot_b, th_b)},
                                  "Le groupe de Galois se lit dans les nombres premiers",
                                  "type de cycles de f modulo p, observé sur les premiers p ≤ 3000 contre les proportions du groupe"),
            "16_quintique_chebotarev")


CAS_AVANCES = {
    "quotients": ("Centre, quotients et homomorphismes", cas_quotients),
    "classification": ("Les 24 groupes d'ordre ≤ 12", cas_classification),
    "sylow": ("Sylow : compter les sous-groupes", cas_sylow),
    "cayley": ("Cayley : tout groupe est un groupe de permutations", cas_cayley),
    "burnside": ("Burnside : colorier à symétrie près", cas_burnside),
    "crypto": ("Diffie–Hellman et RSA dans (Z/nZ)*", cas_crypto),
    "rubik": ("Le Rubik's Cube : 43 trillions d'éléments", cas_rubik),
    "polynomes": ("Le groupe de Galois d'un polynôme (degrés 2-4)", cas_polynomes),
    "quintique": ("Prouver qu'un quintique n'est pas résoluble", cas_quintique),
}
