"""
Globale Assemblierung fuer Q4-Elemente.

Dasselbe Prinzip wie stab/assembly.py: Knoten i belegt die globalen DOFs
[2*i, 2*i+1]. Ein Element hat jetzt 4 Knoten (n0, n1, n2, n3) statt 2 und
damit 8 lokale DOFs:

    [2*n0, 2*n0+1, 2*n1, 2*n1+1, 2*n2, 2*n2+1, 2*n3, 2*n3+1]

Die Verschiebungen eines Elements holst du aus dem globalen Vektor u und
bringst sie in die Form (4, 2), die q4/element.py erwartet:

    u_e = u[dofs].reshape(4, 2)
"""
import numpy as np

from q4 import element


def dofs_of_element(conn) -> list:
    """Globale DOF-Indizes eines Elements mit den Knoten conn = (n0, n1, n2, n3)."""
    raise NotImplementedError("TODO: dofs_of_element")


def assemble(nodes: np.ndarray, elements: list, u: np.ndarray,
             lam: float, mu: float, t: float):
    """Baut den globalen internen Kraftvektor F_int und die globale Tangente K_T.

    nodes:    (n_nodes, 2)  Referenzkoordinaten
    elements: Liste von (n0, n1, n2, n3), gegen den Uhrzeigersinn
    u:        (2*n_nodes,)  aktueller globaler Verschiebungsvektor
    lam, mu:  Lame-Parameter (fuer alle Elemente gleich)
    t:        Dicke der Scheibe

    Rueckgabe: F_int (2*n_nodes,), K_T (2*n_nodes, 2*n_nodes)
    """
    raise NotImplementedError("TODO: assemble (Q4)")
