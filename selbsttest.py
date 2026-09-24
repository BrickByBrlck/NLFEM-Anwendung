"""
Selbst-Check fuer alle Module. Nach jeder Code-Aenderung laufen lassen.

Kein Ersatz fuer main.py am Ende, sondern ein schneller Weg, JEDE Funktion
einzeln zu pruefen, waehrend du sie schreibst -- ohne auf main.py zu warten,
das erst ganz am Ende laeuft (und vorher an jedem fehlenden Stueck abbricht).

Ausfuehren: py selbsttest.py

Gibt zu jedem Check OK/FEHLER aus. Ein FEHLER sagt dir, WAS nicht passt
(Wert oder Form/Typ), nicht die richtige Loesung -- die findest du selbst.
Noch nicht implementierte Funktionen werden als "noch offen" markiert und
ueberspringen den Rest ihres Blocks, damit ein spaeter Fehler nicht die
frueheren Checks verdeckt.
"""
import numpy as np

import material
import element
import assembly
import newton_raphson
import load_stepping


def check(name, condition, hint=""):
    status = "OK   " if condition else "FEHLER"
    print(f"[{status}] {name}")
    if not condition and hint:
        print(f"         -> {hint}")


def try_block(title, func):
    print(f"\n--- {title} ---")
    try:
        func()
    except NotImplementedError:
        print("  (noch offen -- ueberspringe)")
    except Exception as e:
        print(f"  [FEHLER] Ausnahme beim Ausfuehren: {e!r}")


def test_material():
    S = material.stress(E_tt=0.01, Emod=100.0)
    check("stress(0.01, 100) ist eine einzelne Zahl (kein Array)",
          np.ndim(S) == 0,
          "stress() soll eine Zahl zurueckgeben, keine Liste/Array.")
    check("stress(0.01, 100) == 1.0",
          abs(S - 1.0) < 1e-9,
          "Formel ist S_tt = Emod * E_tt.")

    Ct1 = material.tangent_modulus(E_tt=0.0, Emod=100.0)
    Ct2 = material.tangent_modulus(E_tt=5.0, Emod=100.0)
    check("tangent_modulus ist fuer verschiedene E_tt konstant (lineares Material)",
          abs(Ct1 - Ct2) < 1e-12,
          "dS_tt/dE_tt haengt bei S_tt = Emod*E_tt nicht von E_tt ab.")
    check("tangent_modulus(_, 100) == 100",
          abs(Ct1 - 100.0) < 1e-9)


def test_element_strain_undeformed():
    X1, X2 = np.array([0.0, 0.0]), np.array([2.0, 0.0])
    u1, u2 = np.array([0.0, 0.0]), np.array([0.0, 0.0])
    E = element.bar_strain(X1, X2, u1, u2)
    check("bar_strain gibt bei u=0 eine einzelne Zahl zurueck (kein Array)",
          np.ndim(E) == 0,
          "d ist ein 2er-Vektor. d.d (Skalarprodukt, np.dot(d,d) bzw. d@d) ist "
          "eine Zahl -- d**2 ist elementweises Quadrieren und bleibt ein "
          "2er-Array. Das ist ein sehr haeufiger NumPy-Fehler an dieser Stelle.")
    if np.ndim(E) == 0:
        check("bar_strain(unverformter Stab) == 0.0",
              abs(E - 0.0) < 1e-12,
              "Ohne Verschiebung ist d gleich e, also d.d == 1, also E_tt == 0.")


def test_element_strain_small_axial_stretch():
    # Stab entlang der x-Achse, Knoten 2 wird ein kleines Stueck WEITER
    # entlang des Stabs verschoben -> reine Dehnung, kein Querversatz.
    L0 = 2.0
    delta = 1e-3
    X1, X2 = np.array([0.0, 0.0]), np.array([L0, 0.0])
    u1 = np.array([0.0, 0.0])
    u2 = np.array([delta, 0.0])
    E = element.bar_strain(X1, X2, u1, u2)
    if np.ndim(E) != 0:
        print("  (ueberspringe Zahlenvergleich, Form stimmt noch nicht -- siehe oben)")
        return
    E_linear_approx = delta / L0
    check("bar_strain bei kleiner reiner Laengsdehnung ~ delta/L0 (linearisierter Grenzfall)",
          abs(E - E_linear_approx) < 1e-4,
          f"erwartet ungefaehr {E_linear_approx:.6f}, bekommen {E}.")


def test_tangent_vs_finite_differences():
    # Der wichtigste Einzeltest in der nichtlinearen FEM: K_e muss die
    # Ableitung von f_int nach u sein. Wir pruefen das numerisch mit
    # zentralen Differenzen in einem beliebigen, deutlich verformten Zustand.
    X1, X2 = np.array([0.0, 0.0]), np.array([1.0, 0.5])
    u = np.array([0.03, -0.02, -0.05, -0.2])
    Emod, A, h = 1000.0, 1.0, 1e-7

    def f(v):
        return element.internal_force(X1, X2, v[:2], v[2:], Emod, A)

    K = element.tangent_stiffness(X1, X2, u[:2], u[2:], Emod, A)
    check("tangent_stiffness hat Form (4, 4)", np.shape(K) == (4, 4))
    K_fd = np.zeros((4, 4))
    for j in range(4):
        du = np.zeros(4)
        du[j] = h
        K_fd[:, j] = (f(u + du) - f(u - du)) / (2 * h)
    err = np.abs(K - K_fd).max()
    check("tangent_stiffness == finite Differenzen von internal_force",
          err < 1e-4,
          f"max. Abweichung {err:.3e}. Meist fehlt/stimmt K_mat oder K_geo nicht.")
    check("tangent_stiffness ist symmetrisch", np.allclose(K, K.T))


# Referenzbeispiel = das Von-Mises-Fachwerk aus main.py
NODES = np.array([[0.0, 0.0], [1.0, 0.5], [2.0, 0.0]])
ELEMENTS = [(0, 1), (1, 2)]
FREE_DOFS = [2, 3]


def test_assembly():
    F_int, K_T = assembly.assemble(NODES, ELEMENTS, np.zeros(6), 1000.0, 1.0)
    check("assemble gibt (F_int, K_T) mit Formen (6,) und (6, 6) zurueck",
          np.shape(F_int) == (6,) and np.shape(K_T) == (6, 6))
    check("F_int ist bei u = 0 gleich null", np.allclose(F_int, 0.0))
    check("K_T ist symmetrisch", np.allclose(K_T, K_T.T))


def test_newton_zero_load():
    u, res = newton_raphson.solve_load_step(NODES, ELEMENTS, np.zeros(6), np.zeros(6),
                                            FREE_DOFS, 1000.0, 1.0)
    check("Nulllast: u bleibt 0", np.allclose(u, 0.0))
    check("Nulllast: genau 1 Newton-Iteration", len(res) == 1,
          f"bekommen: {len(res)} Iterationen")


def test_newton_quadratic_convergence():
    F_ext = np.array([0, 0, 0, -15.0, 0, 0])
    u, res = newton_raphson.solve_load_step(NODES, ELEMENTS, np.zeros(6), F_ext,
                                            FREE_DOFS, 1000.0, 1.0)
    check("Newton konvergiert (||R|| < 1e-10)", res[-1] < 1e-10,
          f"letztes ||R|| = {res[-1]:.3e} nach {len(res)} Iterationen")
    # Konvergenzordnung aus den letzten drei Residuen: ~1 linear, ~2 quadratisch
    r1, r2, r3 = res[-3:]
    order = np.log(r3 / r2) / np.log(r2 / r1)
    check("Konvergenzordnung ist quadratisch (> 1.5)", order > 1.5,
          f"geschaetzte Ordnung {order:.2f}. Tangente passt nicht zu f_int.")


def test_reference_result():
    # Bekanntes Ergebnis des Von-Mises-Fachwerks bei F = -15, 20 Lastschritte.
    # Aendert sich diese Zahl, hat eine Code-Aenderung das Ergebnis veraendert.
    F_max = np.array([0, 0, 0, -15.0, 0, 0])
    loads, disps, _ = load_stepping.run(NODES, ELEMENTS, F_max, FREE_DOFS,
                                        1000.0, 1.0, n_steps=20)
    check("load_stepping gibt Formen (21,) und (21, 6) zurueck",
          np.shape(loads) == (21,) and np.shape(disps) == (21, 6))
    uy = disps[-1, 3]
    check("Referenzwert: uy(Knoten 1) = -0.048853 bei F = -15",
          abs(uy - (-0.048853)) < 1e-6,
          f"bekommen: {uy:.6f}")


if __name__ == "__main__":
    try_block("material.py", test_material)
    try_block("element.py: unverformter Stab", test_element_strain_undeformed)
    try_block("element.py: kleine reine Laengsdehnung", test_element_strain_small_axial_stretch)
    try_block("element.py: Tangente vs. finite Differenzen", test_tangent_vs_finite_differences)
    try_block("assembly.py", test_assembly)
    try_block("newton_raphson.py: Nulllast", test_newton_zero_load)
    try_block("newton_raphson.py: quadratische Konvergenz", test_newton_quadratic_convergence)
    try_block("load_stepping.py: Referenzergebnis", test_reference_result)
    print("\nFertig. Alles OK? Dann committen. Sonst erst den FEHLER beheben.")
