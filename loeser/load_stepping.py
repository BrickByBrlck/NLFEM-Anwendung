"""
Inkrementelle Laststeigerung (load stepping) um `solve_load_step` herum.

Warum ueberhaupt in Schritten laden, statt F_max in einem Schritt aufzubringen?
Newton-Raphson braucht eine Startschaetzung nahe der Loesung, um zu
konvergieren. Bei kleinen Lastschritten ist die Loesung des vorigen Schritts
eine gute Startschaetzung fuer den naechsten -- das ist der ganze Trick.
"""
import numpy as np

from loeser.newton_raphson import solve_load_step


def run(nodes, elements, F_max, free_dofs, Emod, A, n_steps=20):
    """Steigert die Last linear von 0 auf F_max in n_steps Schritten.

    F_max: (2*n_nodes,) globaler Lastvektor bei voller Last (die meisten
           Eintraege sind 0, nur der/die belasteten DOFs sind ungleich 0)

    Rueckgabe:
      load_history:  (n_steps+1,) angelegte Kraft an `track_dof` je Schritt
                      (track_dof = np.argmax(np.abs(F_max)), also automatisch
                      der am staerksten belastete Freiheitsgrad)
      disp_history:  (n_steps+1, 2*n_nodes) Verschiebung nach jedem Lastschritt
      last_residuals: Liste der ||R_frei||-Werte des LETZTEN Lastschritts
                       (fuer den Konvergenz-Plot)
    """
    track_dof=np.argmax(np.abs(F_max))
    u=np.zeros(2*len(nodes))
    load_history = [0.0]
    disp_history = [u.copy()]

    for step in range(1,n_steps+1):
        F_ext = (step / n_steps) * F_max
        u, res_norms= solve_load_step(nodes,elements,u,F_ext,free_dofs,Emod,A)

        load_history.append(F_ext[track_dof])
        disp_history.append(u.copy())

    return np.array(load_history), np.array(disp_history), res_norms


# ---------------------------------------------------------------------------
# Projekt 3: dieselbe Laststeigerung, aber ohne Wissen ueber das Element
# ---------------------------------------------------------------------------

def run_general(assemble_fn, n_dofs, F_max, free_dofs, n_steps=20, tol=1e-8):
    """Wie run(), aber mit assemble_fn statt (nodes, elements, Emod, A).

    Ruft pro Lastschritt newton_raphson.solve(assemble_fn, u, F_ext, free_dofs,
    tol=tol) auf (solve oben in der Import-Zeile ergaenzen). assemble_fn: siehe dort.

    n_dofs: Laenge des globalen Verschiebungsvektors (frueher 2*len(nodes))

    Rueckgabe:
      disp_history:   (n_steps+1, n_dofs) Verschiebung nach jedem Lastschritt,
                       Zeile 0 ist der unbelastete Zustand
      last_residuals: Liste der ||R_frei||-Werte des LETZTEN Lastschritts
    """
    raise NotImplementedError("TODO: run_general")
