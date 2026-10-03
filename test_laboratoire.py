"""Tests des cas avancés : python -m pytest -q"""

import importlib.util
import math
import pathlib
import random

import pytest

import arithmetique as ar
import catalogue as cat
import cas_avances
import denombrement as dn
import galois_polynomes as gp
import permgroup as pg
import quiz
from groupes import Groupe


# --- Groupe : centre, quotient, produit, homomorphismes, Sylow ---------------
def test_centre_et_equation_des_classes():
    D4, Q8, S3 = cat.diedral(4), cat.quaternions(), cat.symetrique(3)
    assert len(D4.centre()) == 2 and len(Q8.centre()) == 2 and len(S3.centre()) == 1
    for G in (D4, Q8, S3, cat.symetrique(4), cat.alterne(5)):
        z, classes = G.equation_des_classes()
        assert z + sum(classes) == len(G)
        assert all(len(G) % c == 0 for c in classes)
    assert cat.cyclique(6).centre() == list(cat.cyclique(6).elements)


def test_centralisateur_et_normalisateur():
    D4 = cat.diedral(4)
    assert len(D4.centralisateur(D4.elements[1])) == 4          # rotation de 90° : le sous-groupe des rotations
    rotations = D4.engendre(D4.elements[1])
    assert len(D4.normalisateur(rotations)) == 8                # distingué : normalisateur = tout
    S4 = cat.symetrique(4)
    H = S4.engendre(cat.depuis_cycles(4, (1, 2)))
    assert len(S4.normalisateur(H)) == 4


def test_quotients():
    D4, V4 = cat.diedral(4), cat.klein()
    assert D4.quotient(D4.centre()).isomorphisme(V4)
    S4 = cat.symetrique(4)
    klein = [cat.depuis_cycles(4)] + [cat.depuis_cycles(4, *c) for c in
                                      (((1, 2), (3, 4)), ((1, 3), (2, 4)), ((1, 4), (2, 3)))]
    assert S4.quotient(klein).isomorphisme(cat.symetrique(3))
    assert cat.alterne(4).quotient(klein).isomorphisme(cat.cyclique(3))
    with pytest.raises(ValueError):
        S4.quotient(S4.engendre(cat.depuis_cycles(4, (1, 2))))   # pas distingué


def test_produit_direct():
    P = cat.cyclique(2).produit_direct(cat.cyclique(3))
    assert P.verifier().est_groupe and len(P) == 6 and P.est_cyclique()
    assert P.isomorphisme(cat.cyclique(6))
    assert not cat.cyclique(2).produit_direct(cat.cyclique(2)).est_cyclique()


def test_homomorphismes():
    S4, Z2, Z12 = cat.symetrique(4), cat.cyclique(2), cat.cyclique(12)
    signe = lambda p: 0 if cat.signature_permutation(p) == 1 else 1
    assert S4.est_homomorphisme(signe, Z2)
    ker = S4.noyau(signe, Z2)
    assert len(ker) == 12 and len(S4) == len(ker) * len(S4.image(signe))
    assert Z12.est_homomorphisme(lambda x: x % 4, cat.cyclique(4))
    assert not Z12.est_homomorphisme(lambda x: x * x % 12, Z12)
    assert Z12.est_homomorphisme({a: 0 for a in Z12.elements}, Z12)          # morphisme trivial (dict)


def test_groupes_simples():
    assert cat.alterne(5).est_simple() and cat.cyclique(7).est_simple()
    for G in (cat.cyclique(6), cat.alterne(4), cat.symetrique(4), cat.diedral(5), cat.quaternions()):
        assert not G.est_simple()
    assert not cat.cyclique(1).est_simple()


@pytest.mark.parametrize("G, p, taille, nb", [
    (cat.symetrique(4), 2, 8, 3), (cat.symetrique(4), 3, 3, 4), (cat.alterne(4), 2, 4, 1),
    (cat.alterne(5), 5, 5, 6), (cat.symetrique(5), 2, 8, 15), (cat.symetrique(5), 3, 3, 10),
])
def test_sylow(G, p, taille, nb):
    P, n = G.sylow(p)
    assert len(P) == taille and n == nb
    assert G.sous_groupe(P).verifier().est_groupe
    assert n % p == 1 and (len(G) // taille) % n == 0


def test_converse_de_lagrange_fausse():
    assert 6 not in {len(H) for H in cat.alterne(4).sous_groupes()}


def test_cayley():
    for G in (cat.cyclique(5), cat.symetrique(3), cat.diedral(4), cat.quaternions()):
        rep = G.representation_reguliere()
        assert len(set(rep.values())) == len(G)
        assert all(rep[G.op(a, b)] == cat.composer(rep[a], rep[b]) for a in G for b in G)
        assert pg.GroupePermutations(list(rep.values())).ordre() == len(G)


# --- Catalogue : nouveaux groupes -----------------------------------------------
def test_groupes_d_ordre():
    attendu = [1, 1, 1, 2, 1, 2, 1, 5, 2, 2, 1, 5]
    for n in range(1, 13):
        groupes = cat.groupes_d_ordre(n)
        assert len(groupes) == attendu[n - 1]
        for i, G in enumerate(groupes):
            assert len(G) == n and G.verifier().est_groupe
            assert all(G.isomorphisme(H) is None for H in groupes[i + 1:])
    with pytest.raises(ValueError):
        cat.groupes_d_ordre(13)


def test_nouveaux_groupes():
    assert cat.dicyclique(2).isomorphisme(cat.quaternions())
    D3 = cat.dicyclique(3)
    assert len(D3) == 12 and D3.verifier().est_groupe and not D3.est_abelien()
    assert D3.isomorphisme(cat.diedral(6)) is None and D3.isomorphisme(cat.alterne(4)) is None
    assert cat.gl2_mod(2).isomorphisme(cat.symetrique(3))
    assert len(cat.sl2_mod(3)) == 24 and cat.sl2_mod(3).verifier().est_groupe
    assert cat.sl2_mod(3).isomorphisme(cat.symetrique(4)) is None
    assert len(cat.gl2_mod(3)) == 48
    A = cat.groupe_affine(5)
    assert len(A) == 20 and A.verifier().est_groupe and A.est_resoluble()
    assert len(set(cat.permutations_du_groupe_affine(5))) == 20


# --- Groupes de permutations et Rubik ------------------------------------------------
def test_ordre_schreier_sims_petits_groupes():
    S5 = pg.GroupePermutations([cat.depuis_cycles(5, (1, 2)), cat.depuis_cycles(5, (1, 2, 3, 4, 5))])
    assert S5.ordre() == 120 == len(S5.enumerer())
    A5 = pg.GroupePermutations([cat.depuis_cycles(5, (1, 2, 3)), cat.depuis_cycles(5, (1, 2, 3, 4, 5))])
    assert A5.ordre() == 60 == len(A5.enumerer())
    assert A5.contient(cat.depuis_cycles(5, (1, 3, 5))) and not A5.contient(cat.depuis_cycles(5, (1, 2)))
    assert S5.contient(cat.depuis_cycles(5, (1, 2)))
    D4 = pg.GroupePermutations([cat.rotation(4, 1), cat.reflexion(4, 1)])
    assert D4.ordre() == 8 and D4.est_transitif()
    assert pg.GroupePermutations([cat.depuis_cycles(6, (1, 2, 3), (4, 5))]).ordre() == 6
    assert pg.GroupePermutations([tuple(range(4))]).ordre() == 1


def test_generateurs_invalides():
    with pytest.raises(ValueError):
        pg.GroupePermutations([])
    with pytest.raises(ValueError):
        pg.GroupePermutations([(0, 0, 1)])


def test_outils_permutations():
    p = cat.depuis_cycles(7, (1, 2, 3), (4, 5))
    assert pg.ordre_permutation(p) == 6 and pg.type_de_cycles(p) == (3, 2, 1, 1)
    assert pg.compose(p, pg.inverse(p)) == pg.identite(7)
    assert pg.compose(p, p) == cat.composer(p, p)


def test_rubik_ordres_connus():
    assert pg.cube_rubik(2).ordre() == 3_674_160 == math.factorial(7) * 3 ** 6
    cube = pg.cube_rubik(3)
    assert cube.ordre() == 43_252_003_274_489_856_000
    assert cube.ordre() == math.factorial(8) * 3 ** 8 * math.factorial(12) * 2 ** 12 // 12
    assert ar.facteurs_premiers(cube.ordre()) == {2: 27, 3: 14, 5: 3, 7: 2, 11: 1}
    assert len(cube.base) >= 10


def test_rubik_mouvements():
    R, U = pg.mouvement_cube("R"), pg.mouvement_cube("U")
    assert pg.ordre_permutation(R) == 4 and pg.ordre_permutation(pg.compose(R, U)) == 105
    assert pg.GroupePermutations([R, U]).ordre() == 73_483_200
    assert pg.ordre_permutation(pg.compose(pg.compose(R, U), pg.compose(pg.inverse(R), pg.inverse(U)))) == 6
    for f in "RLUDFB":
        m = pg.mouvement_cube(f)
        assert sorted(m) == list(range(54)) and pg.ordre_permutation(m) == 4
    assert len(pg.autocollants(3)) == 54 and len(pg.autocollants(2)) == 24


def test_rubik_appartenance_et_aleatoire():
    cube = pg.cube_rubik(3)
    R, U = pg.mouvement_cube("R"), pg.mouvement_cube("U")
    assert cube.contient(pg.compose(R, U))
    echange_deux_coins = list(range(54))
    echange_deux_coins[0], echange_deux_coins[1] = 1, 0                 # permute 2 autocollants : impossible
    assert not cube.contient(echange_deux_coins)
    for _ in range(5):
        assert cube.contient(cube.element_aleatoire())
    d = cube.distribution_des_ordres(100)
    assert sum(d.values()) == 100 and max(d) <= 1260


def test_rotations_du_cube():
    G, produit = pg.rotations_cube()
    assert len(G) == 24
    groupe = dn.groupe_des_rotations_du_cube()
    assert groupe.verifier().est_groupe and groupe.isomorphisme(cat.symetrique(4))
    for genre, n in (("faces", 6), ("sommets", 8), ("aretes", 12), ("diagonales", 4)):
        assert len(pg.points_du_cube(genre)) == n
        assert len({pg.action_sur(genre, m) for m in G}) == 24      # action fidèle
    with pytest.raises(ValueError):
        pg.points_du_cube("tetraedre")


# --- Arithmétique et cryptographie -------------------------------------------------------
def test_phi_et_ordres():
    assert [ar.phi(n) for n in (1, 2, 9, 10, 3233)] == [1, 1, 6, 4, 3120]
    for n in range(2, 60):
        assert ar.phi(n) == len(cat.inversibles_mod(n))
    assert ar.ordre_multiplicatif(5, 23) == 22
    with pytest.raises(ValueError):
        ar.ordre_multiplicatif(4, 6)


def test_racines_primitives():
    assert ar.racines_primitives(23)[:3] == [5, 7, 10]
    assert len(ar.racines_primitives(23)) == ar.phi(22)
    assert ar.plus_petite_racine_primitive(1_000_003) == 2
    for p in (7, 11, 13, 101):
        assert ar.plus_petite_racine_primitive(p) == ar.racines_primitives(p)[0]
    with pytest.raises(ValueError):
        ar.racines_primitives(15)


def test_cyclicite_de_zn_etoile():
    for n in range(2, 61):
        assert cat.inversibles_mod(n).est_cyclique() == ar.a_une_racine_primitive(n), n


def test_euler_et_carmichael():
    for n in (15, 21, 3233, 100):
        f, lam = ar.phi(n), ar.carmichael(n)
        assert f % lam == 0
        for a in range(2, 50):
            if math.gcd(a, n) == 1:
                assert pow(a, f, n) == 1 and pow(a, lam, n) == 1
    assert ar.carmichael(8) == 2 and ar.carmichael(15) == 4


def test_diffie_hellman():
    rng = random.Random(1)
    for _ in range(20):
        p = rng.choice([11, 13, 17, 19, 23, 101])
        g = ar.plus_petite_racine_primitive(p)
        a, b = rng.randint(2, p - 2), rng.randint(2, p - 2)
        A, B, ka, kb = ar.echange_diffie_hellman(p, g, a, b)
        assert ka == kb == pow(g, a * b, p) and A == pow(g, a, p)


def test_logarithme_discret():
    p, g = 1009, ar.plus_petite_racine_primitive(1009)
    for x in (0, 1, 17, 500, 1007):
        h = pow(g, x, p)
        assert ar.logarithme_discret(g, h, p)[0] == x % (p - 1)
        assert pow(g, ar.logarithme_naif(g, h, p)[0], p) == h
    assert ar.logarithme_discret(2, 3, 7)[0] is None                    # 3 n'est pas une puissance de 2 mod 7
    assert ar.logarithme_discret(g, 5, p)[1] < ar.logarithme_naif(g, pow(g, 1000, p), p)[1]


def test_rsa():
    n, f, d, c, m = ar.rsa_jouet(61, 53, 17, 65)
    assert (n, f, d, c, m) == (3233, 3120, 2753, 2790, 65)
    for msg in range(2, 200):
        assert ar.rsa_jouet(61, 53, 17, msg)[4] == msg
    with pytest.raises(ValueError):
        ar.rsa_jouet(61, 53, 15, 65)                                    # e non premier avec φ(n)


def test_premiers():
    assert ar.est_premier(1_000_003) and not ar.est_premier(1) and not ar.est_premier(91)
    assert ar.facteurs_premiers(360) == {2: 3, 3: 2, 5: 1}
    assert ar.est_premier(ar.nombre_premier_aleatoire(24))


# --- Burnside ---------------------------------------------------------------------------------
@pytest.mark.parametrize("n, k, colliers, bracelets", [(4, 2, 6, 6), (5, 2, 8, 8), (6, 2, 14, 13),
                                                         (6, 3, 130, 92), (7, 2, 20, 18)])
def test_colliers_et_bracelets(n, k, colliers, bracelets):
    assert dn.colliers(n, k) == colliers and dn.bracelets(n, k) == bracelets
    assert dn.orbites_force_brute([cat.rotation(n, j) for j in range(n)], n, k) == colliers
    assert dn.orbites_force_brute(list(cat.diedral(n).elements), n, k) == bracelets


def test_coloriages_du_cube():
    assert [dn.coloriages_du_cube("faces", k) for k in (1, 2, 3, 4)] == [1, 10, 57, 240]
    assert dn.coloriages_du_cube("sommets", 2) == 23 and dn.coloriages_du_cube("aretes", 2) == 218
    perms = [pg.action_sur("faces", m) for m in pg.rotations_cube()[0]]
    assert dn.orbites_force_brute(perms, 6, 3) == 57


def test_indice_des_cycles():
    perms = [pg.action_sur("faces", m) for m in pg.rotations_cube()[0]]
    assert dn.indice_des_cycles(perms) == {(4, 1, 1): 6, (3, 3): 8, (2, 2, 2): 6, (2, 2, 1, 1): 3,
                                           (1, 1, 1, 1, 1, 1): 1}


# --- Groupes de Galois de polynômes ---------------------------------------------------------------
@pytest.mark.parametrize("f, d", [([1, 0, -2], 8), ([1, 0, 0, -2], -108), ([1, 0, -3, 1], 81),
                                  ([1, 0, -1, -1], -23), ([1, 0, 0, 0, -2], -2048),
                                  ([1, 1, 1, 1, 1], 125), ([1, 0, 0, -1, -1], -283),
                                  ([1, 0, 0, 0, -1, -1], 2869), ([1, 0, 0, 0, 0, -2], 50000)])
def test_discriminants(f, d):
    assert gp.discriminant(f) == d


@pytest.mark.parametrize("f, nom, ordre", [
    ([1, 0, -2], "Z/2Z", 2), ([1, 0, -1], "trivial", 1),
    ([1, 0, 0, -2], "S3", 6), ([1, 0, -1, -1], "S3", 6), ([1, 0, -3, 1], "Z/3Z", 3),
    ([1, 0, 0, -1], "Z/2Z", 2), ([1, -1, -2, 2], "Z/2Z", 2), ([1, -6, 11, -6], "trivial", 1),
    ([1, 0, 0, 0, -2], "D4", 8), ([1, 0, 0, 0, 1], "V4", 4), ([1, 1, 1, 1, 1], "Z/4Z", 4),
    ([1, 0, 0, 8, 12], "A4", 12), ([1, 0, 0, -1, -1], "S4", 24), ([1, 0, -5, 0, 6], "V4", 4),
    ([1, 0, -4, 0, 2], "Z/4Z", 4), ([1, 0, 4, 0, 2], "Z/4Z", 4), ([1, 0, -5, 0, 5], "Z/4Z", 4),
    ([1, 0, 0, 0, -1], "Z/2Z", 2), ([1, 0, -1, 0, -2], "V4", 4),
])
def test_groupe_de_galois(f, nom, ordre):
    r = gp.groupe_de_galois(f)
    assert (r["nom"], r["ordre"]) == (nom, ordre) and r["raisons"]
    assert len(gp.groupe_catalogue(nom)) == ordre


def test_groupe_de_galois_erreurs():
    for f in ([2, 0, -1], [1, 0, 0, 0, 0, -1], [1, 2, 1], [1, -3, 3, -1]):
        with pytest.raises(ValueError):
            gp.groupe_de_galois(f)


def test_racines_rationnelles():
    assert gp.racines_rationnelles([1, -1, -2, 2]) == [1]
    assert gp.racines_rationnelles([2, -3, 1]) == [gp.Fraction(1, 2), gp.Fraction(1)]
    assert gp.racines_rationnelles([1, 0, 0, 0]) == [0]
    assert gp.racines_rationnelles([1, 0, 1]) == []


def test_types_de_cycles_mod_p():
    f = [1, 0, 0, 0, -1, -1]
    assert gp.types_de_cycles_mod_p(f, 3) == (5,) and gp.types_de_cycles_mod_p(f, 2) == (3, 2)
    assert gp.types_de_cycles_mod_p(f, 19) is None and gp.types_de_cycles_mod_p(f, 151) is None
    # x³ − 2 : p ≡ 2 (mod 3) ⇒ x ↦ x³ est bijective, une seule racine ; p ≡ 1 (mod 3) ⇒ 0 ou 3 racines
    assert gp.types_de_cycles_mod_p([1, 0, 0, -2], 5) == (2, 1)
    assert gp.types_de_cycles_mod_p([1, 0, 0, -2], 7) == (3,)            # 2 n'est pas un cube modulo 7
    assert gp.types_de_cycles_mod_p([1, 0, 0, -2], 31) == (1, 1, 1)      # 2^10 ≡ 1 : 2 est un cube modulo 31
    for p in gp.premiers_jusqu_a(200):
        t = gp.types_de_cycles_mod_p(f, p)
        assert t is None or sum(t) == 5


def test_theoreme_de_dedekind_sur_un_groupe_connu():
    # x^3 − 2 : groupe S3 ; types possibles 1+1+1, 2+1, 3
    types = set(gp.statistique_de_frobenius([1, 0, 0, -2], 400)[0])
    assert types <= {(1, 1, 1), (2, 1), (3,)} and {(2, 1), (3,)} <= types
    # x^2 + 1 : p ≡ 1 (mod 4) ⇒ deux racines, p ≡ 3 (mod 4) ⇒ irréductible
    for p in (5, 13, 17):
        assert gp.types_de_cycles_mod_p([1, 0, 1], p) == (1, 1)
    for p in (3, 7, 11):
        assert gp.types_de_cycles_mod_p([1, 0, 1], p) == (2,)


def test_preuve_s5():
    ok, p5, ptr = gp.prouve_s5([1, 0, 0, 0, -1, -1])
    assert ok and p5 == 3 and ptr == 163
    ok2, _, ptr2 = gp.prouve_s5([1, 0, 0, 0, 0, -2], borne=1500)
    assert not ok2 and ptr2 is None                                      # x^5 − 2 : jamais de transposition


def test_chebotarev():
    S5 = list(cat.symetrique(5))
    th = gp.frequences_theoriques(S5)
    assert th[(5,)] == pytest.approx(24 / 120) and th[(2, 1, 1, 1)] == pytest.approx(10 / 120)
    obs, total = gp.statistique_de_frobenius([1, 0, 0, 0, -1, -1], 2500)
    for t in ((5,), (4, 1), (3, 2), (3, 1, 1)):
        assert abs(obs.get(t, 0) / total - th[t]) < 0.06
    th20 = gp.frequences_theoriques(cat.permutations_du_groupe_affine(5))
    obs20, total20 = gp.statistique_de_frobenius([1, 0, 0, 0, 0, -2], 2500)
    assert set(obs20) <= set(th20) and abs(obs20[(4, 1)] / total20 - 0.5) < 0.08


# --- Quiz avancé ---------------------------------------------------------------------------------------
def test_questions_avancees_coherentes():
    rng = random.Random(7)
    for f in quiz.QUESTIONS[9:]:
        for _ in range(12):
            q = f(rng)
            assert q.enonce and q.reponse and q.explication and q.indice
    assert len(quiz.QUESTIONS) == 17


def test_quiz_rsa_et_sylow():
    q = quiz.q_rsa(random.Random(0))
    assert int(q.reponse) in (40, 60, 72, 120, 64)
    q = quiz.q_sylow(random.Random(1))
    n = int(q.reponse)
    assert n >= 1 and "Sylow" in q.explication


# --- Cas avancés et programme principal --------------------------------------------------------------------
@pytest.mark.parametrize("cle", list(cas_avances.CAS_AVANCES))
def test_chaque_cas_avance(cle, capsys):
    import matplotlib
    matplotlib.use("Agg")
    figures = []

    def montrer(fabrique, nom):
        figures.append((nom, fabrique()))

    cas_avances.CAS_AVANCES[cle][1](montrer)
    sortie = capsys.readouterr().out
    assert sortie.strip() and "Traceback" not in sortie
    assert "✘" not in sortie.replace("pas un générateur", "")
    import matplotlib.pyplot as plt
    for nom, fig in figures:
        assert fig.get_axes(), nom
        plt.close(fig)


def _programme():
    chemin = pathlib.Path(__file__).with_name("groups-theory-EGalois.py")
    spec = importlib.util.spec_from_file_location("programme_lab", chemin)
    programme = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(programme)
    return programme


def test_programme_a_seize_cas(capsys):
    programme = _programme()
    assert len(programme.CAS) == 16
    programme.main(["--cas", "rubik", "--cas", "quintique", "--sans-graphique"])
    sortie = capsys.readouterr().out
    assert "43 252 003 274 489 856 000" in sortie and "Gal(f) = S5 : True" in sortie


def test_programme_sauve_les_figures(tmp_path, capsys):
    programme = _programme()
    programme.main(["--cas", "classification", "--cas", "burnside", "--sauver", str(tmp_path)])
    assert (tmp_path / "09_classification.png").exists() and (tmp_path / "12_burnside_cube.png").exists()
