"""
Netz (mesh) fuer eine rechteckige Scheibe aus Q4-Elementen.

Das Rechteck [0, Lx] x [0, Ly] wird in nx Elemente in x-Richtung und ny
Elemente in y-Richtung geteilt. Das ergibt (nx+1) * (ny+1) Knoten.

Knotennummerierung zeilenweise von unten nach oben, in jeder Zeile von links
nach rechts. Knoten in Spalte i (0..nx) und Zeile j (0..ny):

    Nummer = j * (nx + 1) + i          Koordinaten = [i * Lx/nx, j * Ly/ny]

Beispiel nx = 3, ny = 2:

    8 --- 9 ---10 ---11
    |  3  |  4  |  5  |
    4 --- 5 --- 6 --- 7
    |  0  |  1  |  2  |
    0 --- 1 --- 2 --- 3

Jedes Element bekommt seine 4 Knoten GEGEN den Uhrzeigersinn, beginnend unten
links (passend zum Referenzelement in shape.py). Element 0 oben ist also
(0, 1, 5, 4), Element 4 ist (5, 6, 10, 9).

Im Uhrzeigersinn nummerierte Elemente haben det(J) < 0 -- das waere ein
"umgestuelptes" Element und rechnet Unsinn.
"""
import numpy as np


def rect_mesh(nx: int, ny: int, Lx: float, Ly: float):
    """Erzeugt Knoten und Elemente der Rechteckscheibe.

    Rueckgabe:
      nodes:    np.array der Form ((nx+1)*(ny+1), 2) mit [X, Y] je Knoten
      elements: Liste von (n0, n1, n2, n3)-Tupeln, nx*ny Stueck,
                elementweise zeilenweise wie im Bild oben
    """
    raise NotImplementedError("TODO: rect_mesh")
