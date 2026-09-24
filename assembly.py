"""
Globale Assemblierung (siehe `assemble_stab` in U3).

Jeder Knoten hat 2 Freiheitsgrade [ux, uy]. Knoten i belegt die globalen
DOFs [2*i, 2*i+1].

Ein Element verbindet zwei Knoten (n1, n2) -> lokale DOFs
[2*n1, 2*n1+1, 2*n2, 2*n2+1] in genau dieser Reihenfolge (passend zur
Reihenfolge [u1x,u1y,u2x,u2y] in element.py).
"""
import numpy as np

import element


def dofs_of_element(n1: int, n2: int) -> list:
    """Globale DOF-Indizes eines Elements, Reihenfolge [n1x,n1y,n2x,n2y]."""
    return [2 * n1, 2 * n1 + 1, 2 * n2, 2 * n2 + 1]


def assemble(nodes: np.ndarray, elements: list, u: np.ndarray, Emod: float, A: float):
    """Baut den globalen internen Kraftvektor F_int und die globale Tangente K_T.

    nodes:    (n_nodes, 2)  Referenzkoordinaten
    elements: Liste von (n1, n2) Knotenindex-Paaren
    u:        (2*n_nodes,)  aktueller globaler Verschiebungsvektor
    Emod, A:  Materialparameter (hier fuer alle Elemente gleich)

    Rueckgabe: F_int (2*n_nodes,), K_T (2*n_nodes, 2*n_nodes)
    """
    n_nodes=len(nodes)
    n_dofs=2*n_nodes

    F_int=np.zeros(n_dofs)
    K_T= np.zeros((n_dofs,n_dofs))

    for n1,n2 in elements:
        X1= nodes[n1]
        X2= nodes[n2]
        dofs=dofs_of_element(n1,n2)
        u1=u[dofs[0:2]]
        u2=u[dofs[2:4]]
        f_e=element.internal_force(X1, X2, u1, u2, Emod, A)
        k_e=element.tangent_stiffness(X1, X2, u1, u2, Emod, A)
        F_int[dofs]+=f_e
        K_T[np.ix_(dofs,dofs)]+=k_e

    return F_int, K_T
