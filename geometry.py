r"""
Geometrie des 2D-Eiffelturms: Knoten, Staebe, Auflager.

Der Turm besteht aus n_levels Etagen. Jede Etage i (i = 0 ... n_levels) hat
einen linken und einen rechten Knoten auf der Hoehe

    y_i = height * i / n_levels

Ganz oben sitzt zusaetzlich eine Spitze (apex) bei (0, height + spire).

Knotennummerierung (wichtig, darauf bauen elements und base_dofs auf):
    Etage i, links  -> Knoten 2*i
    Etage i, rechts -> Knoten 2*i + 1
    Spitze          -> Knoten 2*(n_levels + 1)   (= letzter Knoten)

Beispiel n_levels = 2:

              6            Spitze
             / \
            4---5          Etage 2
            |\ /|
            | X |
            |/ \|
           2-----3         Etage 1
           |\   /|
           | \ / |
           |  X  |
           | / \ |
           |/   \|
          0       1        Etage 0 = Fusspunkte (fest eingespannt)

Staebe pro Etagenfeld (zwischen Etage i und i+1), mit
L = 2*i, R = 2*i+1, Lo = 2*i+2, Ro = 2*i+3:
    (L, Lo), (R, Ro)      die zwei Beine
    (Lo, Ro)              waagerechter Stab oben
    (L, Ro), (R, Lo)      die zwei Diagonalen (X-Verband)
Dazu oben zwei Staebe von den obersten Knoten zur Spitze.

Warum die Diagonalen? Ein Viereck aus 4 Staeben mit Gelenken ist
beweglich (ein "Mechanismus") -- es kann sich zum Parallelogramm
verschieben, ohne dass sich ein Stab dehnt. Dann ist K_T singulaer und
np.linalg.solve bricht ab. Erst die Diagonalen machen das Fachwerk steif.

Breite: Der echte Eiffelturm wird nach oben ungefaehr exponentiell
schmaler (Eiffel hat die Form auf Windlasten hin optimiert):

    w(y) = base_half_width * exp(-k * y)
    k    = ln(base_half_width / top_half_width) / height

Damit ist w(0) = base_half_width und w(height) = top_half_width.
"""
import numpy as np


def half_width(y, height, base_half_width, top_half_width):
    """Halbe Turmbreite w(y) auf Hoehe y (Formel siehe oben).

    TODO:
      1. k = np.log(base_half_width / top_half_width) / height
      2. return base_half_width * np.exp(-k * y)
    """
    raise NotImplementedError("TODO: half_width")


def eiffel_2d(n_levels=8, height=3.0, base_half_width=0.625, top_half_width=0.05,
              spire=0.3):
    """Erzeugt Knoten und Staebe des 2D-Eiffelturms.

    Rueckgabe:
      nodes:    np.array der Form (2*(n_levels+1) + 1, 2) mit [x, y] je Knoten
      elements: Liste von (n1, n2)-Paaren, 5*n_levels + 2 Stueck

    TODO:
      1. Knoten: leere Liste nodes = []
         fuer i in range(n_levels + 1):
           - y = height * i / n_levels
           - w = half_width(y, height, base_half_width, top_half_width)
           - links [-w, y] und rechts [w, y] anhaengen (genau in der Reihenfolge!)
         zum Schluss die Spitze [0.0, height + spire] anhaengen
         und nodes = np.array(nodes)
      2. Staebe: leere Liste elements = []
         fuer i in range(n_levels):
           - L, R, Lo, Ro = 2*i, 2*i+1, 2*i+2, 2*i+3
           - die 5 Staebe des Etagenfelds anhaengen (siehe Liste oben)
         dann die 2 Staebe zur Spitze: (2*n_levels, apex), (2*n_levels+1, apex)
         mit apex = len(nodes) - 1
      3. return nodes, elements
    """
    raise NotImplementedError("TODO: eiffel_2d")


def base_dofs():
    """Fixierte DOFs: beide Fusspunkte (Knoten 0 und 1) in x und y.

    TODO: Knoten i belegt die DOFs [2*i, 2*i+1] (wie in assembly.py).
          Welche 4 DOFs sind das fuer Knoten 0 und 1? Als Liste zurueckgeben.
    """
    raise NotImplementedError("TODO: base_dofs")
