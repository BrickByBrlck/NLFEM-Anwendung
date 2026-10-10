"""
Ansatzfunktionen (shape functions) und Gauss-Punkte des 4-Knoten-Elements (Q4).

Das Element wird NICHT in seinen echten Koordinaten beschrieben, sondern auf
einem immer gleichen Referenzelement: dem Quadrat [-1, 1] x [-1, 1] mit den
Koordinaten (xi, eta). Die vier Knoten sitzen in den Ecken, gegen den
Uhrzeigersinn nummeriert:

        eta
         ^
    3 ---|--- 2          Knoten 0: (xi, eta) = (-1, -1)
    |    |    |          Knoten 1:             (+1, -1)
    |    +----|--> xi    Knoten 2:             (+1, +1)
    |         |          Knoten 3:             (-1, +1)
    0 ------- 1

Jeder Knoten I hat eine Ansatzfunktion N_I, die an "ihrem" Knoten 1 ist und an
den drei anderen 0 (bilinear dazwischen):

    N_I(xi, eta) = 1/4 * (1 + xi_I * xi) * (1 + eta_I * eta)

xi_I, eta_I sind die Eckkoordinaten von Knoten I (stehen unten in NODE_XI).

Ableitungen nach den Referenzkoordinaten (Produktregel, je ein Faktor faellt weg):

    dN_I/dxi  = 1/4 * xi_I  * (1 + eta_I * eta)
    dN_I/deta = 1/4 * eta_I * (1 + xi_I  * xi)

Gauss-Quadratur 2x2: Ein Integral ueber das Referenzelement wird durch eine
gewichtete Summe an 4 Punkten ersetzt,

    Integral f dxi deta  ~  Summe_g  w_g * f(xi_g, eta_g)

mit den Punkten (+-1/sqrt(3), +-1/sqrt(3)) und allen Gewichten w_g = 1.
"""
import numpy as np

# Eckkoordinaten (xi_I, eta_I) der vier Knoten, Zeile I = Knoten I
NODE_XI = np.array([
    [-1.0, -1.0],
    [+1.0, -1.0],
    [+1.0, +1.0],
    [-1.0, +1.0],
])


def shape_functions(xi: float, eta: float) -> np.ndarray:
    """N (4,): Werte der vier Ansatzfunktionen am Punkt (xi, eta)."""
    xi_I = NODE_XI[:, 0]
    eta_I = NODE_XI[:, 1]

    N=0.25 * (1.0 + xi_I * xi) * (1.0 + eta_I * eta)

    return N


    


def shape_derivatives(xi: float, eta: float) -> np.ndarray:
    """dN_dxi (4, 2): Zeile I = [dN_I/dxi, dN_I/deta] am Punkt (xi, eta)."""

    xi_I = NODE_XI[:, 0]
    eta_I = NODE_XI[:, 1]

    dN_dxi = 0.25 * xi_I * (1 + eta_I * eta)
    dN_deta = 0.25 * eta_I * (1 + xi_I * xi)

    return np.column_stack((dN_dxi,dN_deta))


def gauss_points():
    """2x2-Gauss-Quadratur auf dem Referenzelement.

    Rueckgabe: (points, weights)
      points:  (4, 2)  Zeile g = [xi_g, eta_g]
      weights: (4,)    Gewichte w_g
    Die Reihenfolge der vier Punkte ist egal.
    """

    a = 1.0/np.sqrt(3.0)
    points = np.array([
        [-a, -a],
        [+a, -a],
        [+a, +a],
        [-a, +a]
        ])
    weights = np.array([1.0,1.0,1.0,1.0])
    return points, weights
    