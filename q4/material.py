"""
Materialgesetz fuer das 2D-Kontinuumselement: St. Venant-Kirchhoff (SVK),
ebener Verzerrungszustand (plane strain).

Das ist dasselbe Gesetz wie in U2, nur in 2D aufgeschrieben:

    S = lam * tr(E) * I + 2 * mu * E

Ebener Verzerrungszustand heisst: in Dickenrichtung dehnt sich nichts,
E_33 = 0. Deshalb ist tr(E) = E_11 + E_22, und du kannst komplett mit
2x2-Matrizen rechnen. (S_33 ist nicht null, wird aber fuer das Element nicht
gebraucht.)

Lame-Parameter aus E-Modul und Querkontraktionszahl:

    lam = Emod * nu / ((1 + nu) * (1 - 2*nu))
    mu  = Emod / (2 * (1 + nu))

Voigt-Notation: Eine symmetrische 2x2-Matrix hat nur 3 verschiedene Eintraege.
Man schreibt sie als Vektor, in diesem Projekt immer in der Reihenfolge

    S_voigt = [S_11, S_22, S_12]
    E_voigt = [E_11, E_22, 2*E_12]      <- Achtung: Faktor 2 bei der Scherung

Der Faktor 2 sorgt dafuer, dass S_voigt . E_voigt = S : E bleibt (im
Doppelskalarprodukt kommt E_12 zweimal vor: als E_12 und als E_21).

Die Materialtangente C = dS/dE (ein Tensor 4. Stufe) wird damit zu einer
3x3-Matrix, die S_voigt = C_voigt @ E_voigt erfuellt:

    C_voigt = [[lam + 2*mu,  lam,         0 ],
               [lam,         lam + 2*mu,  0 ],
               [0,           0,           mu]]

Bei SVK ist C konstant -- sie haengt nicht von E ab. Bei Neo-Hooke (spaeter)
aendert sich genau das.
"""
import numpy as np


def lame(Emod: float, nu: float):
    """Rueckgabe (lam, mu) aus E-Modul und Querkontraktionszahl."""
    raise NotImplementedError("TODO: lame")


def stress(E: np.ndarray, lam: float, mu: float) -> np.ndarray:
    """S (2, 2) aus E (2, 2):  S = lam*tr(E)*I + 2*mu*E."""
    raise NotImplementedError("TODO: stress")


def to_voigt(S: np.ndarray) -> np.ndarray:
    """[S_11, S_22, S_12] (3,) aus einer symmetrischen 2x2-Matrix.

    Nur fuer SPANNUNGEN gedacht (ohne Faktor 2).
    """
    raise NotImplementedError("TODO: to_voigt")


def tangent_voigt(lam: float, mu: float) -> np.ndarray:
    """C_voigt (3, 3), Formel siehe oben."""
    raise NotImplementedError("TODO: tangent_voigt")
