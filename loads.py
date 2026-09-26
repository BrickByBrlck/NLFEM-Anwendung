"""
Lasten fuer den Eiffelturm: Wind und Eigengewicht.

Beide Funktionen geben einen globalen Lastvektor F der Laenge 2*n_nodes
zurueck, in derselben Reihenfolge wie u: [u0x, u0y, u1x, u1y, ...].
Also: x-Kraft an Knoten i steht in F[2*i], y-Kraft in F[2*i+1].
Die Vektoren kann man einfach addieren: F = wind + gewicht.

Wind:
    Wind ist in der Hoehe staerker als am Boden. Ueblich ist ein
    Potenzgesetz fuer das Windprofil:

        Gewicht_i = (y_i / H) ** alpha        (alpha ~ 0.16, H = hoechster Knoten)

    Die Gesamtkraft total_force wird im Verhaeltnis dieser Gewichte auf alle
    Knoten verteilt, horizontal in +x (Wind von links):

        F_x,i = total_force * Gewicht_i / Summe(Gewichte)

    Knoten mit y = 0 (Fusspunkte) bekommen keinen Wind.

Eigengewicht:
    Ein Stab wiegt G = rho_g * A * L0 (Gewicht pro Volumen mal Volumen).
    Jeder der zwei Knoten des Stabs traegt die Haelfte, nach unten (-y).
"""
import numpy as np

import element


def wind_load(nodes, total_force, alpha=0.16):
    """Horizontale Windlast, nach oben staerker werdend (Formel siehe oben).

    TODO:
      1. y = nodes[:, 1] und H = y.max()
      2. Gewichte fuer alle Knoten: (y / H) ** alpha
         (Fusspunkte haben y = 0 und damit automatisch Gewicht 0)
      3. F = np.zeros(2 * len(nodes))
      4. F[0::2] = total_force * gewichte / gewichte.sum()
         (F[0::2] waehlt jeden zweiten Eintrag ab 0 aus = alle x-Komponenten)
      5. return F
    """
    raise NotImplementedError("TODO: wind_load")


def self_weight(nodes, elements, rho_g, A):
    """Eigengewicht der Staebe, je zur Haelfte auf die zwei Knoten (-y).

    TODO:
      1. F = np.zeros(2 * len(nodes))
      2. fuer jeden Stab (n1, n2) in elements:
           - L0 = element.reference_length(nodes[n1], nodes[n2])
           - G = rho_g * A * L0
           - an beiden Knoten G/2 von der y-Komponente abziehen
             (welcher Index in F ist die y-Komponente von Knoten n1?)
      3. return F
    """
    raise NotImplementedError("TODO: self_weight")
