"""
Materialgesetz (constitutive law) des Stabelements.

Für den Stab aus U3 ist das Materialgesetz linear: die gesamte geometrische
Nichtlinearitaet steckt in der Verzerrung E_tt (siehe element.py), nicht im
Materialgesetz selbst.

    S_tt = Emod * E_tt

S_tt: 2. Piola-Kirchhoff-Spannung in Stabrichtung
Emod: Elastizitaetsmodul
E_tt: Green-Lagrange-Verzerrung in Stabrichtung
"""


def stress(E_tt: float, Emod: float) -> float:
    """S_tt = Emod * E_tt.

    TODO: implementieren.
    """
    return E_tt*Emod
    raise NotImplementedError("TODO: S_tt = Emod * E_tt")


def tangent_modulus(E_tt: float, Emod: float) -> float:
    """dS_tt/dE_tt.

    Fuer das lineare Materialgesetz oben ist das eine Konstante (unabhaengig
    von E_tt) — der Parameter E_tt ist trotzdem Teil der Signatur, weil sich
    das bei einem nichtlinearen Materialgesetz (z.B. spaeter Plastizitaet)
    aendern wuerde.

    TODO: implementieren.
    """
    return Emod
