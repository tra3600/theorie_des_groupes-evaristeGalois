"""Tests : python -m pytest -q"""

import pytest

import catalogue as cat
import quiz
from groupes import Groupe


@pytest.mark.parametrize("G", [
    cat.cyclique(7), cat.horloge(), cat.inversibles_mod(15), cat.symetrique(4),
    cat.alterne(4), cat.diedral(5), cat.klein(), cat.quaternions(),
])
def test_vrais_groupes(G):
    assert G.verifier().est_groupe
    assert G.is_group()  # ancienne API conservée
    for a in G:
        assert G.op(a, G.inverse(a)) == G.neutre
        assert len(G) % G.ordre_element(a) == 0  # Lagrange


@pytest.mark.parametrize("G, axiome", [
    (cat.faux_addition_bornee(), "Fermeture"),
    (cat.faux_soustraction(4), "Associativité"),
    (cat.faux_maximum(), "Inverses"),
    (cat.faux_multiplication_avec_zero(5), "Inverses"),
    (cat.faux_multiplication_mod(6), "Fermeture"),
])
def test_faux_groupes_et_axiome_viole(G, axiome):
    rapport = G.verifier()
    assert not rapport.est_groupe
    assert not next(a for a in rapport.axiomes if a.nom == axiome).valide


def test_multiplication_modulo_premier():
    for n in range(2, 14):
        est_premier = all(n % d for d in range(2, n))
        assert cat.faux_multiplication_mod(n).is_group() == est_premier


def test_ancien_exemple_permutations_de_123():
    # Le programme d'origine : permutations de (1, 2, 3), p1∘p2 indexé à partir de 1.
    import itertools
    G = Groupe(list(itertools.permutations([1, 2, 3])),
               lambda p1, p2: tuple(p1[i - 1] for i in p2))
    assert G.is_group()
    assert G.identity == (1, 2, 3)
    assert len(G.inverses) == 6


def test_cardinaux():
    assert [len(cat.symetrique(n)) for n in range(1, 6)] == [1, 2, 6, 24, 120]
    assert len(cat.alterne(5)) == 60
    assert len(cat.diedral(6)) == 12
    assert len(cat.inversibles_mod(12)) == 4


def test_commutativite():
    assert cat.cyclique(10).est_abelien()
    assert cat.klein().est_abelien()
    for G in (cat.symetrique(3), cat.diedral(4), cat.quaternions()):
        a, b = G.contre_exemple_commutativite()
        assert G.op(a, b) != G.op(b, a)


def test_cyclicite_et_isomorphismes():
    assert cat.cyclique(4).est_cyclique() and not cat.klein().est_cyclique()
    assert cat.inversibles_mod(7).est_cyclique()
    assert cat.inversibles_mod(8).isomorphisme(cat.klein())
    assert cat.cyclique(4).isomorphisme(cat.klein()) is None
    assert cat.diedral(3).isomorphisme(cat.symetrique(3))
    assert cat.diedral(4).isomorphisme(cat.quaternions()) is None
    f = cat.inversibles_mod(5).isomorphisme(cat.cyclique(4))
    G, H = cat.inversibles_mod(5), cat.cyclique(4)
    assert all(f[G.op(a, b)] == H.op(f[a], f[b]) for a in G for b in G)


def test_sous_groupes():
    assert len(cat.symetrique(3).sous_groupes()) == 6
    assert len(cat.diedral(4).sous_groupes()) == 10
    assert len(cat.quaternions().sous_groupes()) == 6
    assert len(cat.symetrique(4).sous_groupes()) == 30
    Q8 = cat.quaternions()
    assert all(Q8.est_distingue(H) for H in Q8.sous_groupes())
    S3 = cat.symetrique(3)
    transposition = S3.engendre(cat.depuis_cycles(3, (1, 2)))
    assert not S3.est_distingue(transposition)
    classes = S3.classes_a_gauche(transposition)
    assert len(classes) == 3 and all(len(c) == 2 for c in classes)


def test_classes_de_conjugaison():
    assert sorted(len(c) for c in cat.symetrique(4).classes_de_conjugaison()) == [1, 3, 6, 6, 8]
    assert sorted(len(c) for c in cat.alterne(5).classes_de_conjugaison()) == [1, 12, 12, 15, 20]


def test_galois_resolubilite():
    assert cat.symetrique(3).est_resoluble()
    assert [len(H) for H in cat.symetrique(4).serie_derivee()] == [24, 12, 4, 1]
    S5 = cat.symetrique(5)
    assert not S5.est_resoluble()
    assert [len(H) for H in S5.serie_derivee()] == [120, 60]
    cinq_cycle = cat.depuis_cycles(5, (1, 2, 3, 4, 5))
    assert len(S5.engendre(cinq_cycle, cat.depuis_cycles(5, (1, 2)))) == 120


def test_permutations_utilitaires():
    p = cat.depuis_cycles(5, (1, 2, 3), (4, 5))
    assert cat.notation_cyclique(p) == "(1 2 3)(4 5)"
    assert cat.signature_permutation(p) == -1
    assert cat.notation_cyclique(tuple(range(4))) == "id"


def test_melanges_parfaits():
    ordre = lambda n, i=False: len(cat.groupe_engendre_par_permutation(cat.melange_parfait(n, i)))
    assert ordre(52) == 8
    assert ordre(52, True) == 52
    # Le mélange extérieur garde la première et la dernière carte en place.
    p = cat.melange_parfait(10)
    assert p[0] == 0 and p[-1] == 9
    with pytest.raises(ValueError):
        cat.melange_parfait(7)


def test_puissances_et_horloge():
    G = cat.horloge()
    assert G.op(9, 5) == 2
    assert G.puissance(5, 5) == 1
    assert G.puissance(1, -3) == 9
    assert G.etiquette(0) == "12h"
    assert [a for a in G if G.ordre_element(a) == 12] == [1, 5, 7, 11]


def test_quiz_toutes_bonnes_reponses():
    questions = quiz.generer(18, graine=1)
    reponses = iter(q.reponse for q in questions)
    sorties = []
    score = quiz.jouer(18, graine=1, entree=lambda _: next(reponses), sortie=sorties.append)
    assert score == 18
    assert any("Héritier" in s for s in sorties)


def test_quiz_indice_et_abandon():
    reponses = iter(["?", "mauvaise", "q"])
    score = quiz.jouer(3, graine=2, entree=lambda _: next(reponses), sortie=lambda _: None)
    assert score == 0


def test_programme_principal_sans_graphique(capsys):
    import importlib.util
    import pathlib
    chemin = pathlib.Path(__file__).with_name("groups-theory-EGalois.py")
    spec = importlib.util.spec_from_file_location("programme", chemin)
    programme = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(programme)
    programme.main(["--tout", "--sans-graphique"])
    sortie = capsys.readouterr().out
    assert "NON résoluble" in sortie
    assert "retour après  8 mélanges" in sortie
