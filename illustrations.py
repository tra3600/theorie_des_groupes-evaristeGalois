"""Les dessins : chaque fonction renvoie une figure matplotlib.

Le programme principal décide ensuite s'il faut l'afficher à l'écran ou
l'enregistrer en PNG (option ``--sauver``).
"""

import math

import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyArrowPatch, Polygon, Rectangle

from catalogue import (appliquer_melange, composer, diedral, nature_symetrie,
                       notation_cyclique)

# Palette catégorielle (ordre fixe, validée pour le daltonisme) et encres neutres.
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
ENCRE = "#0b0b0b"
ENCRE_2 = "#52514e"
GRILLE = "#d9d8d3"
FOND = "#fcfcfb"
BON, MAUVAIS = "#008300", "#e34948"

plt.rcParams.update({
    "figure.facecolor": FOND,
    "axes.facecolor": FOND,
    "axes.edgecolor": GRILLE,
    "text.color": ENCRE,
    "axes.labelcolor": ENCRE_2,
    "xtick.color": ENCRE_2,
    "ytick.color": ENCRE_2,
    "font.size": 10,
})


def _finir(fig, titre, sous_titre=None, bas=0.0):
    """Titre + sous-titre placés à distance fixe (en pouces) du haut, puis mise en page."""
    h = fig.get_figheight()
    fig.text(0.5, 1 - 0.32 / h, titre, ha="center", va="center", fontsize=14,
             fontweight="bold", color=ENCRE)
    marge = 0.55
    if sous_titre:
        fig.text(0.5, 1 - 0.62 / h, sous_titre, ha="center", va="center", fontsize=10,
                 color=ENCRE_2)
        marge = 0.85
    fig.tight_layout(rect=(0, bas, 1, 1 - marge / h))


def _couleur_teinte(i, n):
    """Teintes réparties sur le cercle, douces : utilisées quand le texte porte l'identité."""
    import colorsys
    r, g, b = colorsys.hls_to_rgb((i / max(n, 1)) % 1.0, 0.80, 0.55)
    return (r, g, b)


# ----------------------------------------------------------------------
# Table de Cayley
# ----------------------------------------------------------------------
def table_cayley(G):
    """La table de multiplication colorée : chaque ligne et chaque colonne
    contient chaque élément exactement une fois (un carré latin, comme au sudoku)."""
    n = len(G)
    labs = [G.etiquette(e) for e in G.elements]
    taille = min(0.55 * n + 2.5, 13)
    fig, ax = plt.subplots(figsize=(taille, taille * 0.95))
    for i in range(n):
        for j in range(n):
            k = G.table[i][j]
            couleur = _couleur_teinte(k, n) if k is not None else "#ffffff"
            ax.add_patch(Rectangle((j, n - 1 - i), 1, 1, facecolor=couleur,
                                   edgecolor=FOND, linewidth=2))
            if n <= 24:
                ax.text(j + 0.5, n - 0.5 - i, labs[k] if k is not None else "?",
                        ha="center", va="center", fontsize=max(6, 12 - n // 3), color=ENCRE)
    ax.set_xlim(0, n)
    ax.set_ylim(0, n)
    ax.set_xticks([j + 0.5 for j in range(n)])
    ax.set_xticklabels(labs if n <= 24 else [], rotation=90 if n > 8 else 0)
    ax.set_yticks([n - 0.5 - i for i in range(n)])
    ax.set_yticklabels(labs if n <= 24 else [])
    ax.xaxis.tick_top()
    ax.set_xlabel("b (colonne)")
    ax.set_ylabel("a (ligne)")
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_aspect("equal")
    _finir(fig, f"Table de Cayley de {G.nom} : case (a, b) = a·b",
           "Chaque ligne et chaque colonne contient chaque élément une seule fois — un vrai sudoku !", bas=0)
    return fig


# ----------------------------------------------------------------------
# Graphe de Cayley
# ----------------------------------------------------------------------
def graphe_cayley(G, generateurs, disposition="circulaire"):
    """Un point par élément ; une flèche x → x·g pour chaque générateur g.

    Suivre les flèches, c'est « se promener » dans le groupe."""
    graphe = nx.MultiDiGraph()
    graphe.add_nodes_from(G.elements)
    for c, g in enumerate(generateurs):
        for x in G.elements:
            graphe.add_edge(x, G.op(x, g), gen=c)

    if disposition == "circulaire":
        pos = nx.circular_layout(graphe)
    elif disposition == "anneaux":
        # Anneau extérieur : e, g, g², ... (g = premier générateur).  Chaque anneau
        # suivant place x·h au même angle que x (h = dernier générateur) : les
        # flèches de h deviennent des rayons et ne croisent plus rien.
        g, h = generateurs[0], generateurs[-1]
        anneau, x = [], G.neutre
        while x not in anneau:
            anneau.append(x)
            x = G.op(x, g)
        nb_anneaux = len(G) // len(anneau)
        angles = {k: math.pi / 2 - 2 * math.pi * k / len(anneau) for k in range(len(anneau))}
        pos, niveau = {}, 0
        while True:
            r = 1.0 - (0.55 * niveau / (nb_anneaux - 1) if nb_anneaux > 1 else 0)
            for k, x in enumerate(anneau):
                pos.setdefault(x, (r * math.cos(angles[k]), r * math.sin(angles[k])))
            anneau = [G.op(x, h) for x in anneau]
            niveau += 1
            if len(pos) == len(G) or niveau >= nb_anneaux:
                break
        for x in G.elements:  # filet de sécurité si h ne suffit pas à tout atteindre
            pos.setdefault(x, (0.0, 0.0))
    else:
        pos = nx.spring_layout(nx.Graph(graphe), seed=7, iterations=200)

    fig, ax = plt.subplots(figsize=(8, 8))
    for c, g in enumerate(generateurs):
        aretes = [(u, v) for u, v, d in graphe.edges(data=True) if d["gen"] == c and u != v]
        involution = G.ordre_element(g) == 2
        nx.draw_networkx_edges(
            graphe, pos, edgelist=aretes, ax=ax, edge_color=SERIES[c], width=2,
            arrows=True, arrowstyle="-" if involution else "-|>", arrowsize=16, node_size=1300,
            connectionstyle="arc3,rad=0.12" if not involution else "arc3,rad=0.0")
    nx.draw_networkx_nodes(graphe, pos, ax=ax, node_size=1300, node_color=FOND,
                           edgecolors=ENCRE_2, linewidths=1.5)
    nx.draw_networkx_labels(graphe, pos, labels={x: G.etiquette(x) for x in G.elements},
                            ax=ax, font_size=9, font_color=ENCRE)
    poignees = [plt.Line2D([], [], color=SERIES[c], lw=2,
                           label=f"× {G.etiquette(g)}" + ("  (aller-retour)" if G.ordre_element(g) == 2 else ""))
                for c, g in enumerate(generateurs)]
    ax.legend(handles=poignees, loc="lower center", bbox_to_anchor=(0.5, -0.08),
              ncol=len(generateurs), frameon=False)
    ax.set_axis_off()
    ax.set_aspect("equal")
    _finir(fig, f"Graphe de Cayley de {G.nom}",
           "Chaque flèche colorée = multiplier à droite par un générateur", bas=0.02)
    return fig


# ----------------------------------------------------------------------
# Symétries d'un polygone
# ----------------------------------------------------------------------
def _sommet(n, i, r=1.0):
    a = math.pi / 2 + 2 * math.pi * i / n
    return r * math.cos(a), r * math.sin(a)


def _dessiner_polygone(ax, p, titre):
    n = len(p)
    ax.add_patch(Polygon([_sommet(n, i) for i in range(n)], closed=True,
                         facecolor="#eef3fb", edgecolor=ENCRE_2, linewidth=1.5))
    genre, k = nature_symetrie(p)
    if genre == "miroir":
        a = math.pi / 2 + math.pi * k / n
        ax.plot([-1.25 * math.cos(a), 1.25 * math.cos(a)],
                [-1.25 * math.sin(a), 1.25 * math.sin(a)],
                linestyle="--", color=SERIES[1], linewidth=2)
    elif k:
        debut, fin = math.pi / 2 + 0.3, math.pi / 2 + 2 * math.pi * k / n - 0.3
        ts = [debut + (fin - debut) * t / 40 for t in range(41)]
        xs, ys = [0.38 * math.cos(t) for t in ts], [0.38 * math.sin(t) for t in ts]
        ax.plot(xs, ys, color=SERIES[0], linewidth=2)
        ax.annotate("", xy=(xs[-1], ys[-1]), xytext=(xs[-3], ys[-3]),
                    arrowprops=dict(arrowstyle="-|>", color=SERIES[0], lw=2))
    # Le sommet qui était en position i se retrouve en position p[i].
    for i in range(n):
        x, y = _sommet(n, p[i])
        ax.scatter([x], [y], s=380, color=SERIES[i % len(SERIES)], edgecolors=FOND,
                   linewidths=2, zorder=3)
        ax.text(x, y, str(i + 1), ha="center", va="center", color="white",
                fontsize=10, fontweight="bold", zorder=4)
    ax.set_title(titre, fontsize=10, color=ENCRE)
    ax.set_xlim(-1.35, 1.35)
    ax.set_ylim(-1.35, 1.35)
    ax.set_aspect("equal")
    ax.set_axis_off()


def symetries_polygone(n):
    """Les 2n façons de reposer un polygone régulier dans son contour."""
    G = diedral(n)
    fig, axes = plt.subplots(2, n, figsize=(2.3 * n, 5.2))
    for ax, p in zip(axes.flat, G.elements):
        genre, _ = nature_symetrie(p)
        legende = "rotation" if genre == "rotation" else "miroir"
        _dessiner_polygone(ax, p, f"{G.etiquette(p)}  ({legende})\n{notation_cyclique(p)}")
    _finir(fig, f"Les {2 * n} symétries du polygone à {n} côtés : le groupe D{n}",
           "Ligne du haut : rotations (flèche bleue)   ·   Ligne du bas : miroirs (axe orange en pointillés)", bas=0)
    return fig


def composition_polygone(n, a, b):
    """Montre a∘b en trois étapes : départ, après b, puis après a."""
    G = diedral(n)
    ident = G.neutre
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.9))
    _dessiner_polygone(axes[0], ident, "Départ")
    _dessiner_polygone(axes[1], b, f"On applique {G.etiquette(b)}")
    _dessiner_polygone(axes[2], composer(a, b),
                       f"Puis {G.etiquette(a)}  ⇒  même effet que {G.etiquette(G.op(a, b))}")
    _finir(fig, f"{G.etiquette(a)} ∘ {G.etiquette(b)} = {G.etiquette(G.op(a, b))}",
           "Composer deux symétries donne encore une symétrie : c'est la fermeture", bas=0)
    return fig


# ----------------------------------------------------------------------
# Horloge
# ----------------------------------------------------------------------
def horloge(depart, duree):
    """Le cadran : partir de ``depart`` heures et avancer de ``duree`` heures."""
    fig, ax = plt.subplots(figsize=(6, 6.4))
    ax.add_patch(plt.Circle((0, 0), 1.05, facecolor="#f4f3ef", edgecolor=ENCRE_2, linewidth=2))
    for h in range(12):
        a = math.pi / 2 - 2 * math.pi * h / 12
        ax.text(0.88 * math.cos(a), 0.88 * math.sin(a), str(h or 12),
                ha="center", va="center", fontsize=14,
                fontweight="bold" if h in (depart % 12, (depart + duree) % 12) else "normal",
                color=ENCRE)
    arrivee = (depart + duree) % 12
    # Trajet parcouru : une spirale qui rentre d'un cran à chaque tour complet.
    etapes = max(1, abs(duree)) * 12
    xs, ys = [], []
    for t in range(etapes + 1):
        heures = duree * t / etapes
        a = math.pi / 2 - 2 * math.pi * (depart + heures) / 12
        r = 0.66 - 0.08 * abs(heures) / 12
        xs.append(r * math.cos(a))
        ys.append(r * math.sin(a))
    ax.plot(xs, ys, color=SERIES[1], linewidth=3)
    if len(xs) > 1:
        ax.annotate("", xy=(xs[-1], ys[-1]), xytext=(xs[-2], ys[-2]),
                    arrowprops=dict(arrowstyle="-|>", color=SERIES[1], lw=3))
    for h, c in ((depart, SERIES[0]), (arrivee, SERIES[2])):
        a = math.pi / 2 - 2 * math.pi * h / 12
        ax.plot([0, 0.45 * math.cos(a)], [0, 0.45 * math.sin(a)], color=c, linewidth=5,
                solid_capstyle="round")
    ax.scatter([0], [0], s=60, color=ENCRE)
    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(-1.2, 1.2)
    ax.set_aspect("equal")
    ax.set_axis_off()
    signe = "+" if duree >= 0 else "−"
    _finir(fig, f"{depart or 12}h {signe} {abs(duree)}h = {arrivee or 12}h",
           f"Aiguille bleue : départ · orange : trajet · verte : arrivée.  Dans Z/12Z, {depart} {signe} {abs(duree)} ≡ {arrivee}", bas=0)
    return fig


# ----------------------------------------------------------------------
# Mélange de cartes
# ----------------------------------------------------------------------
def melanges(p, titre_melange):
    """Chaque ligne = le paquet après un mélange de plus.  Couleur = numéro de la carte."""
    n = len(p)
    paquet, lignes = list(range(n)), [list(range(n))]
    while True:
        paquet = appliquer_melange(paquet, p)
        lignes.append(paquet)
        if paquet == lignes[0]:
            break
    k = len(lignes) - 1
    rampe = LinearSegmentedColormap.from_list("bleus", ["#dbe9fa", "#2a78d6", "#0d2f5c"])
    hauteur = min(1.8 + 0.28 * len(lignes), 16)
    fig, ax = plt.subplots(figsize=(11, hauteur))
    ax.imshow(lignes, cmap=rampe, aspect="auto", interpolation="nearest")
    ax.set_yticks(range(len(lignes)))
    ax.set_yticklabels(["départ"] + [f"après {i}" for i in range(1, len(lignes))],
                       fontsize=8 if len(lignes) > 20 else 10)
    ax.set_xlabel("position dans le paquet (0 = dessus)")
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)
    _finir(fig, f"{titre_melange} de {n} cartes : retour à l'ordre initial après {k} mélanges",
           "Couleur = numéro de la carte (clair = dessus du paquet d'origine, foncé = dessous)", bas=0)
    return fig


# ----------------------------------------------------------------------
# Treillis des sous-groupes
# ----------------------------------------------------------------------
def treillis_sous_groupes(G):
    """Diagramme de Hasse : un nœud par sous-groupe, rangés par taille.

    Les sous-groupes distingués (normaux) sont en vert."""
    sgs = [frozenset(H) for H in G.sous_groupes()]
    par_ordre = {}
    for H in sgs:
        par_ordre.setdefault(len(H), []).append(H)
    ordres = sorted(par_ordre)
    pos = {}
    for niveau, o in enumerate(ordres):
        ligne = par_ordre[o]
        for i, H in enumerate(ligne):
            pos[H] = (i - (len(ligne) - 1) / 2, niveau)

    fig, ax = plt.subplots(figsize=(max(7, 1.4 * max(len(v) for v in par_ordre.values())), 1.3 * len(ordres) + 2))
    for H in sgs:
        # Arêtes de couverture : H ⊂ K sans intermédiaire
        au_dessus = [K for K in sgs if H < K]
        for K in au_dessus:
            if not any(H < M < K for M in au_dessus):
                ax.plot(*zip(pos[H], pos[K]), color=GRILLE, linewidth=1.5, zorder=1)
    for H in sgs:
        distingue = G.est_distingue(H)
        if len(H) == 1:
            nom = "{e}"
        elif len(H) == len(G):
            nom = G.nom
        else:
            gens = _petits_generateurs(G, H)
            nom = "⟨" + ", ".join(G.etiquette(g) for g in gens) + "⟩"
        ax.text(*pos[H], nom, ha="center", va="center", fontsize=9, zorder=3,
                color=ENCRE,
                bbox=dict(boxstyle="round,pad=0.35", facecolor="#e3f2e3" if distingue else FOND,
                          edgecolor=BON if distingue else ENCRE_2, linewidth=1.5))
    ax.set_yticks(range(len(ordres)))
    ax.set_yticklabels([f"ordre {o}" for o in ordres])
    ax.set_xticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)
    xmax = max(abs(x) for x, _ in pos.values()) + 0.8
    ax.set_xlim(-xmax, xmax)
    ax.set_ylim(-0.6, len(ordres) - 0.4)
    _finir(fig, f"Les {len(sgs)} sous-groupes de {G.nom}",
           "Encadré vert = sous-groupe distingué (normal).  Un trait relie H au-dessous de K quand H ⊂ K.", bas=0)
    return fig


def _petits_generateurs(G, H):
    """Quelques éléments qui engendrent H (glouton)."""
    H = set(H)
    gens, engendre = [], {G.neutre}
    for x in sorted(H, key=lambda x: -G.ordre_element(x)):
        if x not in engendre:
            gens.append(x)
            engendre = set(G.engendre(*gens))
        if engendre == H:
            break
    return gens


# ----------------------------------------------------------------------
# Séries dérivées (Galois)
# ----------------------------------------------------------------------
NOMS_DERIVES = {("S3", 1): "A3", ("S4", 1): "A4", ("S4", 2): "V4", ("S5", 1): "A5"}


def series_derivees(groupes):
    """Pour chaque groupe, la chaîne G ⊇ G' ⊇ G'' ... : descend-elle jusqu'à {e} ?"""
    fig, ax = plt.subplots(figsize=(10, 1.1 * len(groupes) + 1.2))
    for ligne, G in enumerate(groupes):
        serie = G.serie_derivee()
        y = len(groupes) - 1 - ligne
        resoluble = len(serie[-1]) == 1
        for etape, H in enumerate(serie):
            x = etape * 2.2
            ax.add_patch(Rectangle((x, y - 0.3), 1.6, 0.6, facecolor=FOND,
                                   edgecolor=ENCRE_2, linewidth=1.5))
            nom = G.nom if etape == 0 else ("{e}" if len(H) == 1 else
                                            NOMS_DERIVES.get((G.nom, etape), G.nom + "'" * etape))
            ax.text(x + 0.8, y + 0.07, nom, ha="center", va="center", fontsize=11, fontweight="bold")
            ax.text(x + 0.8, y - 0.16, f"{len(H)} élément{'s' if len(H) > 1 else ''}", ha="center", va="center",
                    fontsize=8, color=ENCRE_2)
            if etape:
                ax.annotate("", xy=(x - 0.05, y), xytext=(x - 0.55, y),
                            arrowprops=dict(arrowstyle="-|>", color=ENCRE_2, lw=1.5))
        xf = len(serie) * 2.2
        if resoluble:
            verdict, couleur = "✔ résoluble : formule par radicaux", BON
        else:
            dernier = NOMS_DERIVES.get((G.nom, len(serie) - 1), G.nom + "'" * (len(serie) - 1))
            verdict, couleur = f"✘ bloqué : {dernier} est son propre dérivé", MAUVAIS
        ax.text(xf - 0.2, y, verdict, va="center", fontsize=10, color=couleur, fontweight="bold")
    ax.set_xlim(-0.3, 15)
    ax.set_ylim(-0.7, len(groupes) - 0.3)
    ax.set_axis_off()
    _finir(fig, "Galois : une équation est résoluble par radicaux ⇔ son groupe est résoluble",
           "On prend les commutateurs aba⁻¹b⁻¹ encore et encore : arrive-t-on à {e} ?", bas=0)
    return fig


# ----------------------------------------------------------------------
# Figures du laboratoire avancé
# ----------------------------------------------------------------------
def puissances_modulaires(p, gs):
    """Pour chaque g de ``gs`` : les points 1..p−1 sur un cercle, une flèche de x vers g·x.

    Un générateur trace un seul grand cycle ; sinon la multiplication se découpe en cycles plus courts."""
    gs = [gs] if isinstance(gs, int) else list(gs)
    fig, axes = plt.subplots(1, len(gs), figsize=(5.6 * len(gs), 5.8))
    axes = [axes] if len(gs) == 1 else list(axes)
    pts = {x: (math.sin(2 * math.pi * x / p), math.cos(2 * math.pi * x / p)) for x in range(1, p)}
    for ax, g in zip(axes, gs):
        ordre, x = 1, g % p
        while x != 1:
            x = x * g % p
            ordre += 1
        for x in range(1, p):
            y = x * g % p
            if y == x:
                continue
            ax.add_patch(FancyArrowPatch(pts[x], pts[y], arrowstyle="-|>", mutation_scale=11,
                                         color=SERIES[0], lw=1.3, shrinkA=13, shrinkB=13,
                                         connectionstyle="arc3,rad=0.15"))
        orbite = set()
        x = 1
        while x not in orbite:
            orbite.add(x)
            x = x * g % p
        for x, (px, py) in pts.items():
            ax.scatter([px], [py], s=520, color=SERIES[2] if x in orbite else GRILLE, zorder=3,
                       edgecolor=ENCRE_2)
            ax.text(px, py, str(x), ha="center", va="center", fontsize=10, zorder=4)
        ax.set_xlim(-1.35, 1.35)
        ax.set_ylim(-1.35, 1.35)
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_title(f"× {g} : ordre {ordre}" + (" = p − 1 : générateur ✔" if ordre == p - 1
                                                  else " < p − 1 : pas un générateur"),
                     fontsize=11, color=BON if ordre == p - 1 else MAUVAIS)
    _finir(fig, f"Multiplier modulo {p}", "vert : l'orbite de 1, c'est-à-dire les puissances de g")
    return fig


def barres_burnside(series, titre, sous_titre=None):
    """series : {étiquette: (valeurs par k)} ; barres groupées, k = nombre de couleurs."""
    ks = sorted(next(iter(series.values())))
    fig, ax = plt.subplots(figsize=(9, 5))
    largeur = 0.8 / len(series)
    for i, (nom, vals) in enumerate(series.items()):
        xs = [j + i * largeur for j in range(len(ks))]
        barres = ax.bar(xs, [vals[k] for k in ks], largeur, label=nom, color=SERIES[i % len(SERIES)])
        for b in barres:
            ax.text(b.get_x() + b.get_width() / 2, b.get_height(), f"{int(b.get_height())}",
                    ha="center", va="bottom", fontsize=8, color=ENCRE_2)
    ax.set_xticks([j + 0.4 - largeur / 2 for j in range(len(ks))], [f"{k} couleur{'s' if k > 1 else ''}" for k in ks])
    ax.set_yscale("log")
    ax.set_ylabel("coloriages différents (échelle log)")
    ax.legend(frameon=False)
    ax.grid(axis="y", color=GRILLE)
    ax.spines[["top", "right"]].set_visible(False)
    _finir(fig, titre, sous_titre)
    return fig


def distribution_des_ordres(distribution, titre, sous_titre=None):
    """Histogramme des ordres d'éléments tirés au hasard : {ordre: effectif}."""
    fig, ax = plt.subplots(figsize=(10, 4.8))
    ordres = sorted(distribution)
    total = sum(distribution.values())
    ax.bar(range(len(ordres)), [100 * distribution[o] / total for o in ordres], color=SERIES[0])
    ax.set_xticks(range(len(ordres)), [str(o) for o in ordres], rotation=60, fontsize=8)
    ax.set_xlabel("ordre de l'élément")
    ax.set_ylabel("% des éléments tirés")
    ax.grid(axis="y", color=GRILLE)
    ax.spines[["top", "right"]].set_visible(False)
    _finir(fig, titre, sous_titre)
    return fig


def classification(groupes_par_ordre):
    """Tableau des groupes d'ordre 1 à 12 avec leurs invariants (✔/✘)."""
    lignes, couleurs = [], []
    for n, groupes in groupes_par_ordre.items():
        for G in groupes:
            lignes.append([str(n), G.nom, "oui" if G.est_abelien() else "non",
                           "oui" if G.est_cyclique() else "non", str(len(G.centre())),
                           str(len(G.sous_groupes())), "oui" if G.est_resoluble() else "non"])
            couleurs.append(FOND if n % 2 else "#eeeeea")
    fig, ax = plt.subplots(figsize=(9.5, 0.3 * len(lignes) + 1.4))
    ax.axis("off")
    tab = ax.table(cellText=lignes, colLabels=["ordre", "groupe", "abélien", "cyclique", "|centre|",
                                               "sous-groupes", "résoluble"],
                   cellLoc="center", bbox=[0, 0, 1, 1])
    tab.auto_set_font_size(False)
    tab.set_fontsize(9)
    tab.scale(1, 1.25)
    for (i, j), cell in tab.get_celld().items():
        cell.set_edgecolor(GRILLE)
        if i == 0:
            cell.set_facecolor(SERIES[0])
            cell.get_text().set_color("white")
            cell.get_text().set_fontweight("bold")
        else:
            cell.set_facecolor(couleurs[i - 1])
    _finir(fig, "Les 24 groupes d'ordre ≤ 12", "il n'y en a que 24 : chacun a sa carte d'identité")
    return fig


def frobenius(statistiques, titre, sous_titre=None):
    """statistiques : {nom: (observé {type: n}, total, théorique {type: proportion})}."""
    types = sorted({t for obs, _, th in statistiques.values() for t in list(obs) + list(th)},
                   key=lambda t: (-len(t), t), reverse=True)
    fig, axes = plt.subplots(1, len(statistiques), figsize=(6.4 * len(statistiques), 4.8), sharey=True)
    axes = [axes] if len(statistiques) == 1 else list(axes)
    for ax, (nom, (obs, total, th)) in zip(axes, statistiques.items()):
        xs = range(len(types))
        ax.bar([x - 0.2 for x in xs], [100 * obs.get(t, 0) / total for t in types], 0.4,
               color=SERIES[0], label=f"observé ({total} premiers)")
        ax.bar([x + 0.2 for x in xs], [100 * th.get(t, 0) for t in types], 0.4,
               color=SERIES[1], label="prédit par Chebotarev")
        ax.set_xticks(list(xs), ["+".join(map(str, t)) for t in types], rotation=60, fontsize=8)
        ax.set_title(nom, fontsize=11, fontweight="bold")
        ax.set_xlabel("degrés des facteurs de f modulo p")
        ax.legend(frameon=False, fontsize=8)
        ax.grid(axis="y", color=GRILLE)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("% des nombres premiers")
    _finir(fig, titre, sous_titre)
    return fig
