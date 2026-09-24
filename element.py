"""
Stabelement mit grossen Deformationen (2 Knoten, 2D, 4 Freiheitsgrade).

Referenzkonfiguration: Knotenkoordinaten X1, X2  (je (2,)-Array: [X, Y])
Momentane Verschiebungen:            u1, u2  (je (2,)-Array: [ux, uy])

Kinematik (siehe U3):
    L0   = |X2 - X1|                        Referenzlaenge
    e    = (X2 - X1) / L0                   Einheitsvektor entlang des Stabs (Referenz)
    d    = e + (u2 - u1) / L0                aktuelle Stabrichtung / L0
                                              (d = (x2 - x1) / L0 mit x = X + u;
                                               |d| ist i.A. NICHT 1, das ist Absicht)

Green-Lagrange-Verzerrung in Stabrichtung, exakt (nicht linearisiert):
    E_tt = 0.5 * (d . d - 1.0)
         (aequivalent zu E_tt = (Lcur^2 - L0^2) / (2*L0^2), Lcur = aktuelle Laenge)

Spannung und Stabkraft:
    S_tt = material.stress(E_tt, Emod)
    N    = A * S_tt                          Stabkraft (bezogen auf Referenzflaeche A)

Interner Kraftvektor (4,), Reihenfolge [u1x, u1y, u2x, u2y]:
    f_int = N * [-d, d]        (d als 2er-Vektor: f_int = concat(-N*d, N*d))

Elementtangente (4x4), materieller + geometrischer Anteil:
    Ct    = material.tangent_modulus(E_tt, Emod)
    K_mat = (A * Ct / L0) * outer([-d, d], [-d, d])          # 4x4, [-d,d] als 4er-Vektor
    K_geo = (N / L0) * [[ I, -I],
                         [-I,  I]]                             # I = 2x2-Einheitsmatrix
    K_e   = K_mat + K_geo

Das ist exakt der Stoff aus U3 (E_tt, S_tt, Elementresiduum, k_e = materiell +
geometrisch, `assemble_stab`) — hier baust du ihn selbst, nur mit anderen
Variablennamen/Struktur als in den Kursdateien.
"""
import numpy as np

import material


def reference_length(X1: np.ndarray, X2: np.ndarray) -> float:
    """L0 = |X2 - X1|. Fertig, kein TODO."""
    return float(np.linalg.norm(X2 - X1))


def bar_strain(X1: np.ndarray, X2: np.ndarray, u1: np.ndarray, u2: np.ndarray) -> float:
    """Green-Lagrange-Verzerrung E_tt des Stabs.

    TODO: implementieren (L0, e, d berechnen, dann E_tt = 0.5*(d.d - 1)).
    """
    L0=reference_length(X1,X2)
    e=(X2-X1)/L0
    x1,x2=X1+u1,X2+u2
    d=(x2-x1)/L0
    return 0.5*(d@d -1)

    raise NotImplementedError("TODO: E_tt berechnen")


def internal_force(X1: np.ndarray, X2: np.ndarray, u1: np.ndarray, u2: np.ndarray,
                    Emod: float, A: float) -> np.ndarray:
    """Interner Kraftvektor f_int (4,) in der Reihenfolge [u1x,u1y,u2x,u2y].

    TODO: implementieren (E_tt -> S_tt -> N -> f_int = N*[-d, d]).
    """
    L0=reference_length(X1,X2)
    e=(X2-X1)/L0
    d=e+(u2-u1)/L0
    E_tt=0.5*(d@d -1)
    S_tt=material.stress(E_tt,Emod)
   
    
    N=A*S_tt
    f_int = N*np.concatenate([-d,d])
    return f_int
    raise NotImplementedError("TODO: f_int berechnen")


def tangent_stiffness(X1: np.ndarray, X2: np.ndarray, u1: np.ndarray, u2: np.ndarray,
                       Emod: float, A: float) -> np.ndarray:
    """Elementtangente K_e (4x4) = materieller + geometrischer Anteil.
    Elementtangente (4x4), materieller + geometrischer Anteil:
    Ct    = material.tangent_modulus(E_tt, Emod)
    K_mat = (A * Ct / L0) * outer([-d, d], [-d, d])          # 4x4, [-d,d] als 4er-Vektor
    K_geo = (N / L0) * [[ I, -I],
                         [-I,  I]]                             # I = 2x2-Einheitsmatrix
    K_e   = K_mat + K_geo
    TODO: implementieren (siehe Formeln oben: K_mat + K_geo).
    """
    L0=reference_length(X1,X2)
    e=(X2-X1)/L0
    d=e+(u2-u1)/L0
    E_tt=0.5*(d@d -1)
    S_tt=material.stress(E_tt,Emod)
    N=A*S_tt
    Ct=material.tangent_modulus(E_tt,Emod)
    v = np.concatenate([-d, d])
    K_mat = (A*Ct/L0) * np.outer(v, v)
    I=np.eye(2)
    K_geo = (N/L0)*np.block([[I,-I],[-I,I]])
    return K_mat + K_geo
    raise NotImplementedError("TODO: K_e berechnen")
