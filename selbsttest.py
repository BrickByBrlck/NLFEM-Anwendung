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
import geometry
import loads


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
    load_hist, disps, _ = load_stepping.run(NODES, ELEMENTS, F_max, FREE_DOFS,
                                        1000.0, 1.0, n_steps=20)
    check("load_stepping gibt Formen (21,) und (21, 6) zurueck",
          np.shape(load_hist) == (21,) and np.shape(disps) == (21, 6))
    uy = disps[-1, 3]
    check("Referenzwert: uy(Knoten 1) = -0.048853 bei F = -15",
          abs(uy - (-0.048853)) < 1e-6,
          f"bekommen: {uy:.6f}")


# ---------------------------------------------------------------------------
# Projekt 2: Eiffelturm
# ---------------------------------------------------------------------------

def test_half_width():
    w0 = geometry.half_width(0.0, 3.0, 0.625, 0.05)
    wH = geometry.half_width(3.0, 3.0, 0.625, 0.05)
    wm = geometry.half_width(1.5, 3.0, 0.625, 0.05)
    check("half_width(0) == base_half_width", abs(w0 - 0.625) < 1e-12,
          f"bekommen: {w0}")
    check("half_width(height) == top_half_width", abs(wH - 0.05) < 1e-12,
          f"bekommen: {wH}")
    check("half_width(Mitte) == 0.176777 (exponentiell, nicht linear)",
          abs(wm - 0.176777) < 1e-6,
          f"bekommen: {wm:.6f}. Linear waere 0.3375 -- exp(-k*y) verwendet?")


def test_eiffel_geometry():
    nodes, elements = geometry.eiffel_2d(n_levels=3)
    check("n_levels=3: nodes hat Form (9, 2)", np.shape(nodes) == (9, 2),
          f"bekommen: {np.shape(nodes)}. 2 Knoten pro Etage (4 Etagen) + Spitze.")
    check("n_levels=3: 17 Staebe (5 pro Etagenfeld + 2 zur Spitze)",
          len(elements) == 17, f"bekommen: {len(elements)}")
    if np.shape(nodes) != (9, 2):
        return
    check("Fusspunkte liegen bei y = 0", np.allclose(nodes[:2, 1], 0.0))
    check("Knoten 2*i ist links (x < 0), Knoten 2*i+1 rechts (x > 0)",
          np.all(nodes[0:8:2, 0] < 0) and np.all(nodes[1:8:2, 0] > 0),
          "Reihenfolge pro Etage: erst links, dann rechts anhaengen.")
    check("Turm ist spiegelsymmetrisch", np.allclose(nodes[0:8:2, 0], -nodes[1:8:2, 0]))
    check("Spitze ist der letzte Knoten, bei (0, 3.3)", np.allclose(nodes[-1], [0.0, 3.3]),
          f"bekommen: {nodes[-1]}")
    idx = [n for el in elements for n in el]
    check("alle Stab-Knotennummern liegen zwischen 0 und 8",
          min(idx) >= 0 and max(idx) <= 8)
    check("kein Stab kommt doppelt vor",
          len({tuple(sorted(el)) for el in elements}) == len(elements))


def test_eiffel_stable():
    # Der wichtigste Geometrie-Test: Ist das Fachwerk steif (kein Mechanismus)?
    nodes, elements = geometry.eiffel_2d(n_levels=3)
    fixed = geometry.base_dofs()
    check("base_dofs() == [0, 1, 2, 3]", list(fixed) == [0, 1, 2, 3],
          f"bekommen: {fixed}")
    free = [d for d in range(2 * len(nodes)) if d not in fixed]
    _, K_T = assembly.assemble(nodes, elements, np.zeros(2 * len(nodes)), 1000.0, 1.0)
    K_free = K_T[np.ix_(free, free)]
    rank = np.linalg.matrix_rank(K_free)
    check("K_T (freie DOFs) ist regulaer -> Fachwerk ist steif",
          rank == len(free),
          f"Rang {rank} statt {len(free)}: {len(free) - rank} Bewegungsmoeglichkeit(en) "
          "ohne Stabdehnung. Fehlt eine Diagonale oder ein Querstab?")


def test_loads():
    nodes, elements = geometry.eiffel_2d(n_levels=3)
    F = loads.wind_load(nodes, 10.0)
    check("wind_load hat Laenge 2*n_nodes", np.shape(F) == (18,),
          f"bekommen: {np.shape(F)}")
    if np.shape(F) == (18,):
        check("Wind: Summe der x-Kraefte == total_force", abs(F[0::2].sum() - 10.0) < 1e-9,
              f"bekommen: {F[0::2].sum()}")
        check("Wind: keine y-Kraefte", np.allclose(F[1::2], 0.0))
        check("Wind: Fusspunkte bekommen nichts", np.allclose(F[[0, 2]], 0.0))
        check("Wind: oben staerker als unten", F[2 * 6] > F[2 * 2] > 0,
              "Knoten 6 (Etage 3) muss mehr Wind haben als Knoten 2 (Etage 1).")

    G = loads.self_weight(nodes, elements, rho_g=2.0, A=0.5)
    L_total = sum(element.reference_length(nodes[a], nodes[b]) for a, b in elements)
    check("self_weight: keine x-Kraefte", np.allclose(G[0::2], 0.0))
    check("self_weight: Summe der y-Kraefte == -rho_g*A*(Summe aller L0)",
          abs(G[1::2].sum() + 2.0 * 0.5 * L_total) < 1e-9,
          f"bekommen {G[1::2].sum():.6f}, erwartet {-2.0 * 0.5 * L_total:.6f}")


def test_eiffel_reference():
    # Genau die Einstellungen aus main_eiffel.py
    nodes, elements = geometry.eiffel_2d(n_levels=8)
    fixed = geometry.base_dofs()
    free = [d for d in range(2 * len(nodes)) if d not in fixed]
    F = loads.wind_load(nodes, 10.0) + loads.self_weight(nodes, elements, 1.0, 1.0)
    _, disps, res = load_stepping.run(nodes, elements, F, free, 1000.0, 1.0, n_steps=20)
    ux = disps[-1, -2]
    check("Referenzwert: Spitze ux = 0.469546 (main_eiffel.py-Einstellungen)",
          abs(ux - 0.469546) < 1e-5, f"bekommen: {ux:.6f}")
    check("letzter Lastschritt konvergiert", res[-1] < 1e-10,
          f"letztes ||R|| = {res[-1]:.3e}")


if __name__ == "__main__":
    try_block("material.py", test_material)
    try_block("element.py: unverformter Stab", test_element_strain_undeformed)
    try_block("element.py: kleine reine Laengsdehnung", test_element_strain_small_axial_stretch)
    try_block("element.py: Tangente vs. finite Differenzen", test_tangent_vs_finite_differences)
    try_block("assembly.py", test_assembly)
    try_block("newton_raphson.py: Nulllast", test_newton_zero_load)
    try_block("newton_raphson.py: quadratische Konvergenz", test_newton_quadratic_convergence)
    try_block("load_stepping.py: Referenzergebnis", test_reference_result)
    print("\n========== Projekt 2: Eiffelturm ==========")
    try_block("geometry.py: half_width", test_half_width)
    try_block("geometry.py: eiffel_2d", test_eiffel_geometry)
    try_block("geometry.py: Fachwerk steif? + base_dofs", test_eiffel_stable)
    try_block("loads.py: Wind und Eigengewicht", test_loads)
    try_block("Eiffelturm: Referenzergebnis", test_eiffel_reference)
    print("\nFertig. Alles OK? Dann committen. Sonst erst den FEHLER beheben.")
