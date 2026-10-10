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

from stab import material, element, assembly, geometry, loads
from loeser import newton_raphson, load_stepping
from q4 import shape as shape_q4
from q4 import material as material_2d
from q4 import element as element_q4
from q4 import assembly as assembly_q4
from q4 import mesh


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


# ---------------------------------------------------------------------------
# Projekt 3: Scheibe aus Q4-Elementen
# ---------------------------------------------------------------------------

# Ein absichtlich schiefes Viereck: auf einem sauberen Quadrat heben sich
# manche Fehler in der Jacobi-Matrix gegenseitig auf und bleiben unentdeckt.
XQ = np.array([[0.0, 0.0], [2.0, 0.2], [2.3, 1.5], [-0.2, 1.1]])
UQ = np.array([[0.02, -0.01], [0.05, 0.03], [-0.04, 0.06], [0.01, -0.03]])
UNIT = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]])
LAM, MU = 3.0, 2.0


def rigid_rotation(X, phi):
    """Verschiebungen (n, 2) einer Starrkoerperdrehung um den Ursprung."""
    R = np.array([[np.cos(phi), -np.sin(phi)], [np.sin(phi), np.cos(phi)]])
    return X @ R.T - X


def polygon_area(X):
    x, y = X[:, 0], X[:, 1]
    return 0.5 * (x @ np.roll(y, -1) - y @ np.roll(x, -1))


def count_zero_eigenvalues(K):
    ev = np.linalg.eigvalsh(0.5 * (K + K.T))
    return int(np.sum(np.abs(ev) < 1e-9 * np.abs(ev).max()))


def test_q4_shape_functions():
    N = shape_q4.shape_functions(0.3, -0.7)
    check("shape_functions hat Form (4,)", np.shape(N) == (4,), f"bekommen: {np.shape(N)}")
    check("Summe der vier N ist ueberall 1", abs(np.sum(N) - 1.0) < 1e-12,
          f"bekommen: {np.sum(N)}")
    at_nodes = np.array([shape_q4.shape_functions(xi, eta) for xi, eta in shape_q4.NODE_XI])
    check("N_I ist 1 am eigenen Knoten und 0 an den anderen drei",
          np.allclose(at_nodes, np.eye(4)),
          "Zeile J = N ausgewertet an Knoten J, erwartet die Einheitsmatrix.")

    dN = shape_q4.shape_derivatives(0.3, -0.7)
    check("shape_derivatives hat Form (4, 2)", np.shape(dN) == (4, 2),
          f"bekommen: {np.shape(dN)}")
    check("Ableitungen summieren sich ueber die Knoten zu 0",
          np.allclose(np.sum(dN, axis=0), 0.0),
          "Weil die Summe der N konstant 1 ist, muss die Summe der Ableitungen 0 sein.")
    h = 1e-6
    fd_xi = (shape_q4.shape_functions(0.3 + h, -0.7) - shape_q4.shape_functions(0.3 - h, -0.7)) / (2 * h)
    fd_eta = (shape_q4.shape_functions(0.3, -0.7 + h) - shape_q4.shape_functions(0.3, -0.7 - h)) / (2 * h)
    check("Spalte 0 == Ableitung von shape_functions nach xi", np.allclose(dN[:, 0], fd_xi, atol=1e-8))
    check("Spalte 1 == Ableitung von shape_functions nach eta", np.allclose(dN[:, 1], fd_eta, atol=1e-8),
          "xi und eta vertauscht?")


def test_q4_gauss():
    pts, w = shape_q4.gauss_points()
    check("gauss_points gibt Formen (4, 2) und (4,) zurueck",
          np.shape(pts) == (4, 2) and np.shape(w) == (4,),
          f"bekommen: {np.shape(pts)} und {np.shape(w)}")
    check("Summe der Gewichte == 4 (Flaeche des Referenzquadrats)", abs(np.sum(w) - 4.0) < 1e-12,
          f"bekommen: {np.sum(w)}")
    pts = np.asarray(pts)
    check("Integral von xi^2 ueber das Referenzquadrat == 4/3",
          abs(np.sum(w * pts[:, 0] ** 2) - 4 / 3) < 1e-12,
          f"bekommen: {np.sum(w * pts[:, 0] ** 2):.6f}. Liegen die Punkte bei +-1/sqrt(3)?")
    check("Integral von xi^2 * eta^2 == 4/9",
          abs(np.sum(w * pts[:, 0] ** 2 * pts[:, 1] ** 2) - 4 / 9) < 1e-12)
    check("vier verschiedene Punkte", len({(round(a, 9), round(b, 9)) for a, b in pts}) == 4)


def test_q4_geometry():
    dN_dxi = shape_q4.shape_derivatives(0.3, -0.4)
    J = element_q4.jacobian(UNIT, dN_dxi)
    check("jacobian hat Form (2, 2)", np.shape(J) == (2, 2), f"bekommen: {np.shape(J)}")
    check("Einheitsquadrat: J == 0.5 * Einheitsmatrix", np.allclose(J, 0.5 * np.eye(2)),
          "Das Referenzquadrat hat Kantenlaenge 2, das Einheitsquadrat Kantenlaenge 1.")
    J2 = element_q4.jacobian(np.array([[0.0, 0.0], [4.0, 0.0], [4.0, 1.0], [0.0, 1.0]]), dN_dxi)
    check("Rechteck 4 x 1: J == [[2, 0], [0, 0.5]]", np.allclose(J2, [[2.0, 0.0], [0.0, 0.5]]),
          f"bekommen:\n{J2}\n         Transponiert? J[i, j] = dX_i/dxi_j.")

    dN_dX, detJ = element_q4.physical_gradients(XQ, 0.3, -0.4)
    check("physical_gradients gibt dN_dX mit Form (4, 2) und detJ als Zahl zurueck",
          np.shape(dN_dX) == (4, 2) and np.ndim(detJ) == 0)
    # Ein lineares Feld g(X) = a + b.X muss exakt abgeleitet werden.
    b = np.array([0.7, -1.3])
    g = 2.0 + XQ @ b
    check("Gradient eines linearen Felds wird exakt getroffen",
          np.allclose(dN_dX.T @ g, b),
          f"erwartet {b}, bekommen {dN_dX.T @ g}. dN_dX = dN_dxi @ inv(J)?")
    pts, w = shape_q4.gauss_points()
    area = sum(wg * element_q4.physical_gradients(XQ, xi, eta)[1] for (xi, eta), wg in zip(pts, w))
    check("Summe detJ * w ueber die Gauss-Punkte == Flaeche des Vierecks",
          abs(area - polygon_area(XQ)) < 1e-12,
          f"bekommen {area:.6f}, erwartet {polygon_area(XQ):.6f}")


def test_q4_kinematics():
    dN_dX, _ = element_q4.physical_gradients(XQ, 0.3, -0.4)
    F0 = element_q4.deformation_gradient(np.zeros((4, 2)), dN_dX)
    check("u = 0: F == Einheitsmatrix", np.shape(F0) == (2, 2) and np.allclose(F0, np.eye(2)))
    # Homogene Deformation u = H.X  ->  F = I + H an jedem Punkt
    H = np.array([[0.10, 0.25], [-0.05, 0.20]])
    F = element_q4.deformation_gradient(XQ @ H.T, dN_dX)
    check("homogene Deformation u = H.X: F == I + H", np.allclose(F, np.eye(2) + H),
          f"bekommen:\n{F}\n         Transponiert? F[i, j] = delta_ij + du_i/dX_j.")
    E = element_q4.green_lagrange(np.array([[1.2, 0.3], [0.0, 0.9]]))
    check("green_lagrange([[1.2, 0.3], [0, 0.9]]) == [[0.22, 0.18], [0.18, -0.05]]",
          np.allclose(E, [[0.22, 0.18], [0.18, -0.05]]), f"bekommen:\n{E}")
    F_rot = element_q4.deformation_gradient(rigid_rotation(XQ, 0.7), dN_dX)
    check("Starrkoerperdrehung um 0.7 rad: E == 0",
          np.allclose(element_q4.green_lagrange(F_rot), 0.0, atol=1e-12),
          "Der wichtigste Test auf 'wirklich nichtlinear': mit linearisierter "
          "Verzerrung kaeme hier etwas ungleich 0 heraus.")


def test_material_2d():
    lam, mu = material_2d.lame(1000.0, 0.3)
    check("lame(1000, 0.3) == (576.923077, 384.615385)",
          abs(lam - 576.923077) < 1e-5 and abs(mu - 384.615385) < 1e-5,
          f"bekommen: ({lam:.6f}, {mu:.6f})")
    E = np.array([[0.22, 0.18], [0.18, -0.05]])
    check("stress(E = 0) == 0", np.allclose(material_2d.stress(np.zeros((2, 2)), LAM, MU), 0.0))
    S = material_2d.stress(E, LAM, MU)
    check("stress hat Form (2, 2) und ist symmetrisch",
          np.shape(S) == (2, 2) and np.allclose(S, np.transpose(S)))
    check("stress(E, lam=3, mu=2) == [[1.39, 0.72], [0.72, 0.31]]",
          np.allclose(S, [[1.39, 0.72], [0.72, 0.31]]), f"bekommen:\n{S}")
    Sv = material_2d.to_voigt(S)
    check("to_voigt gibt [S_11, S_22, S_12] zurueck",
          np.shape(Sv) == (3,) and np.allclose(Sv, [S[0, 0], S[1, 1], S[0, 1]]),
          f"bekommen: {Sv}")
    C = material_2d.tangent_voigt(LAM, MU)
    check("tangent_voigt hat Form (3, 3) und ist symmetrisch",
          np.shape(C) == (3, 3) and np.allclose(C, np.transpose(C)))
    E_voigt = np.array([E[0, 0], E[1, 1], 2 * E[0, 1]])
    check("C_voigt @ [E_11, E_22, 2*E_12] == to_voigt(stress(E))",
          np.allclose(C @ E_voigt, Sv),
          f"bekommen {C @ E_voigt}, erwartet {Sv}. Scher-Eintrag mu oder 2*mu?")


def test_q4_b_matrix():
    dN_dX, _ = element_q4.physical_gradients(XQ, 0.3, -0.4)
    F = element_q4.deformation_gradient(UQ, dN_dX)
    B = element_q4.b_matrix(F, dN_dX)
    check("b_matrix hat Form (3, 8)", np.shape(B) == (3, 8), f"bekommen: {np.shape(B)}")

    def E_voigt(u_flat):
        E = element_q4.green_lagrange(element_q4.deformation_gradient(u_flat.reshape(4, 2), dN_dX))
        return np.array([E[0, 0], E[1, 1], 2 * E[0, 1]])

    # B ist die Ableitung von E_voigt nach den 8 Knotenverschiebungen.
    h = 1e-6
    B_fd = np.zeros((3, 8))
    for j in range(8):
        du = np.zeros(8)
        du[j] = h
        B_fd[:, j] = (E_voigt(UQ.reshape(8) + du) - E_voigt(UQ.reshape(8) - du)) / (2 * h)
    err = np.abs(B - B_fd).max()
    row_err = np.abs(B - B_fd).max(axis=1)
    check("b_matrix == Ableitung von [E_11, E_22, 2*E_12] nach u_e",
          err < 1e-7,
          f"max. Abweichung je Zeile: {row_err[0]:.1e}, {row_err[1]:.1e}, {row_err[2]:.1e}. "
          "Die Zeile mit der grossen Zahl ist die falsche.")


def test_q4_element():
    t = 0.5
    f0 = element_q4.internal_force(XQ, np.zeros((4, 2)), LAM, MU, t)
    check("internal_force hat Form (8,)", np.shape(f0) == (8,), f"bekommen: {np.shape(f0)}")
    check("internal_force(u = 0) == 0", np.allclose(f0, 0.0))
    f_rot = element_q4.internal_force(XQ, rigid_rotation(XQ, 0.7), LAM, MU, t)
    check("Starrkoerperdrehung um 0.7 rad: internal_force == 0", np.abs(f_rot).max() < 1e-12,
          f"groesster Eintrag: {np.abs(f_rot).max():.3e}")

    # Absoluter Test: Rechteck 2 x 1, gleichmaessig um 10 % in x gestreckt.
    # F_11 = 1.1, E_11 = 0.105, S_11 = (lam + 2*mu) * E_11, P_11 = F_11 * S_11.
    # Kraft auf die rechte Kante = P_11 * Hoehe * Dicke.
    XR = np.array([[0.0, 0.0], [2.0, 0.0], [2.0, 1.0], [0.0, 1.0]])
    u_stretch = np.column_stack([0.1 * XR[:, 0], np.zeros(4)])
    f = element_q4.internal_force(XR, u_stretch, LAM, MU, t)
    expected = 1.1 * (LAM + 2 * MU) * 0.105 * 1.0 * t
    check(f"10 % Streckung: x-Kraft an der rechten Kante == {expected:.5f}",
          abs(f[2] + f[4] - expected) < 1e-10,
          f"bekommen: {f[2] + f[4]:.5f}. Faktor falsch? Dann detJ, Gewicht oder Dicke pruefen.")

    K = element_q4.tangent_stiffness(XQ, UQ, LAM, MU, t)
    check("tangent_stiffness hat Form (8, 8)", np.shape(K) == (8, 8), f"bekommen: {np.shape(K)}")
    check("tangent_stiffness ist symmetrisch", np.allclose(K, np.transpose(K)))
    h = 1e-6
    K_fd = np.zeros((8, 8))
    for j in range(8):
        du = np.zeros(8)
        du[j] = h
        K_fd[:, j] = (element_q4.internal_force(XQ, (UQ.reshape(8) + du).reshape(4, 2), LAM, MU, t)
                      - element_q4.internal_force(XQ, (UQ.reshape(8) - du).reshape(4, 2), LAM, MU, t)) / (2 * h)
    err = np.abs(K - K_fd).max()
    check("tangent_stiffness == finite Differenzen von internal_force",
          err < 1e-6,
          f"max. Abweichung {err:.3e}. Stimmt der Test bei u = 0, aber hier nicht, "
          "liegt es an K_geo.")
    K0 = element_q4.tangent_stiffness(UNIT, np.zeros((4, 2)), LAM, MU, 1.0)
    check("Einheitsquadrat, u = 0: K[0, 0] == (lam + 2*mu)/3 + mu/3 == 3.0",
          abs(K0[0, 0] - 3.0) < 1e-12, f"bekommen: {K0[0, 0]:.6f}")
    nz = count_zero_eigenvalues(K0)
    check("u = 0: K_e hat genau 3 Null-Eigenwerte (2 Verschiebungen + 1 Drehung)",
          nz == 3, f"bekommen: {nz}")


# Kleines Netz von Hand, damit dieser Test nicht von q4/mesh.py abhaengt.
NODES_Q4 = np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [0.0, 1.0], [1.0, 1.1], [2.0, 1.0]])
ELEMENTS_Q4 = [(0, 1, 4, 3), (1, 2, 5, 4)]


def test_assembly_q4():
    dofs = assembly_q4.dofs_of_element((1, 2, 5, 4))
    check("dofs_of_element((1, 2, 5, 4)) == [2, 3, 4, 5, 10, 11, 8, 9]",
          list(dofs) == [2, 3, 4, 5, 10, 11, 8, 9], f"bekommen: {list(dofs)}")
    F_int, K_T = assembly_q4.assemble(NODES_Q4, ELEMENTS_Q4, np.zeros(12), LAM, MU, 1.0)
    check("assemble gibt (F_int, K_T) mit Formen (12,) und (12, 12) zurueck",
          np.shape(F_int) == (12,) and np.shape(K_T) == (12, 12))
    check("F_int ist bei u = 0 gleich null", np.allclose(F_int, 0.0))
    check("K_T ist symmetrisch", np.allclose(K_T, K_T.T))
    check("Knoten 0 und 2 teilen kein Element: K_T-Block dazwischen ist 0",
          np.allclose(K_T[0:2, 4:6], 0.0))
    k1 = element_q4.tangent_stiffness(NODES_Q4[[0, 1, 4, 3]], np.zeros((4, 2)), LAM, MU, 1.0)
    k2 = element_q4.tangent_stiffness(NODES_Q4[[1, 2, 5, 4]], np.zeros((4, 2)), LAM, MU, 1.0)
    check("Knoten 1 gehoert zu beiden Elementen: K_T[2, 2] ist die Summe beider Beitraege",
          abs(K_T[2, 2] - (k1[2, 2] + k2[0, 0])) < 1e-12,
          "+= statt = beim Einsortieren?")
    u_rot = rigid_rotation(NODES_Q4, 0.7).reshape(12)
    F_rot, _ = assembly_q4.assemble(NODES_Q4, ELEMENTS_Q4, u_rot, LAM, MU, 1.0)
    check("Starrkoerperdrehung des ganzen Netzes: F_int == 0", np.abs(F_rot).max() < 1e-12,
          f"groesster Eintrag: {np.abs(F_rot).max():.3e}. u_e = u[dofs].reshape(4, 2)?")
    nz = count_zero_eigenvalues(K_T)
    check("K_T ohne Randbedingungen hat genau 3 Null-Eigenwerte", nz == 3, f"bekommen: {nz}")


def test_mesh():
    nodes, elements = mesh.rect_mesh(3, 2, 6.0, 2.0)
    check("rect_mesh(3, 2, ...): nodes hat Form (12, 2)", np.shape(nodes) == (12, 2),
          f"bekommen: {np.shape(nodes)}")
    check("rect_mesh(3, 2, ...): 6 Elemente mit je 4 Knoten",
          len(elements) == 6 and all(len(el) == 4 for el in elements))
    if np.shape(nodes) != (12, 2) or len(elements) != 6:
        return
    check("Knoten 6 liegt bei (4, 1)", np.allclose(nodes[6], [4.0, 1.0]), f"bekommen: {nodes[6]}")
    check("Knoten 11 liegt bei (6, 2)", np.allclose(nodes[11], [6.0, 2.0]), f"bekommen: {nodes[11]}")
    check("Element 0 == (0, 1, 5, 4)", tuple(elements[0]) == (0, 1, 5, 4),
          f"bekommen: {tuple(elements[0])}")
    check("Element 4 == (5, 6, 10, 9)", tuple(elements[4]) == (5, 6, 10, 9),
          f"bekommen: {tuple(elements[4])}")
    areas = [polygon_area(nodes[list(el)]) for el in elements]
    check("alle Elemente gegen den Uhrzeigersinn (Flaeche > 0)", min(areas) > 0,
          f"kleinste Flaeche: {min(areas):.3f}")
    check("Summe der Elementflaechen == 12", abs(sum(areas) - 12.0) < 1e-12)
    check("jeder Knoten gehoert zu mindestens einem Element",
          {n for el in elements for n in el} == set(range(12)))


def test_general_solver():
    # Die neuen Funktionen muessen mit Staeben dasselbe liefern wie die alten.
    asm = lambda u: assembly.assemble(NODES, ELEMENTS, u, 1000.0, 1.0)
    F_ext = np.array([0, 0, 0, -15.0, 0, 0])
    u_old, res_old = newton_raphson.solve_load_step(NODES, ELEMENTS, np.zeros(6), F_ext,
                                                    FREE_DOFS, 1000.0, 1.0, tol=1e-10)
    u0 = np.zeros(6)
    u_new, res_new = newton_raphson.solve(asm, u0, F_ext, FREE_DOFS, tol=1e-10)
    check("solve: gleiches u wie solve_load_step (Von-Mises-Fachwerk)",
          np.allclose(u_new, u_old, atol=1e-12))
    check("solve: gleiche Anzahl Iterationen", len(res_new) == len(res_old),
          f"{len(res_new)} statt {len(res_old)}")
    check("solve veraendert das uebergebene u0 nicht", np.allclose(u0, 0.0),
          "u = u0.copy() am Anfang?")

    _, disps_old, _ = load_stepping.run(NODES, ELEMENTS, F_ext, FREE_DOFS, 1000.0, 1.0, n_steps=20)
    disps_new, res = load_stepping.run_general(asm, 6, F_ext, FREE_DOFS, n_steps=20, tol=1e-10)
    check("run_general gibt disp_history mit Form (21, 6) zurueck", np.shape(disps_new) == (21, 6),
          f"bekommen: {np.shape(disps_new)}")
    check("run_general: gleiche Verschiebungen wie run", np.allclose(disps_new, disps_old, atol=1e-10))
    check("run_general: letzter Lastschritt konvergiert", res[-1] < 1e-10)


def test_cantilever_reference():
    # Kragbalken 10 x 1 mit 10 x 2 Elementen, P = 1 nach unten, 5 Lastschritte.
    nodes, elements = mesh.rect_mesh(10, 2, 10.0, 1.0)
    lam, mu = material_2d.lame(1000.0, 0.3)
    n_dofs = 2 * len(nodes)
    left = np.where(np.isclose(nodes[:, 0], 0.0))[0]
    right = np.where(np.isclose(nodes[:, 0], 10.0))[0]
    fixed = [d for n in left for d in (2 * n, 2 * n + 1)]
    free = [d for d in range(n_dofs) if d not in fixed]
    F = np.zeros(n_dofs)
    F[2 * right + 1] = -1.0 / len(right)
    asm = lambda u: assembly_q4.assemble(nodes, elements, u, lam, mu, 1.0)
    disps, res = load_stepping.run_general(asm, n_dofs, F, free, n_steps=5)
    tip = right[len(right) // 2]
    uy = disps[-1, 2 * tip + 1]
    check("letzter Lastschritt konvergiert (||R|| < 1e-8)", res[-1] < 1e-8,
          f"letztes ||R|| = {res[-1]:.3e} nach {len(res)} Iterationen")
    check("hoechstens 6 Newton-Iterationen im letzten Lastschritt", len(res) <= 6,
          f"bekommen: {len(res)}. Mehr deutet auf eine falsche Tangente hin.")
    check("Referenzwert: uy am Balkenende = -2.427587",
          abs(uy - (-2.427587)) < 1e-5, f"bekommen: {uy:.6f}")


if __name__ == "__main__":
    try_block("stab/material.py", test_material)
    try_block("stab/element.py: unverformter Stab", test_element_strain_undeformed)
    try_block("stab/element.py: kleine reine Laengsdehnung", test_element_strain_small_axial_stretch)
    try_block("stab/element.py: Tangente vs. finite Differenzen", test_tangent_vs_finite_differences)
    try_block("stab/assembly.py", test_assembly)
    try_block("loeser/newton_raphson.py: Nulllast", test_newton_zero_load)
    try_block("loeser/newton_raphson.py: quadratische Konvergenz", test_newton_quadratic_convergence)
    try_block("loeser/load_stepping.py: Referenzergebnis", test_reference_result)
    print("\n========== Projekt 2: Eiffelturm ==========")
    try_block("stab/geometry.py: half_width", test_half_width)
    try_block("stab/geometry.py: eiffel_2d", test_eiffel_geometry)
    try_block("stab/geometry.py: Fachwerk steif? + base_dofs", test_eiffel_stable)
    try_block("stab/loads.py: Wind und Eigengewicht", test_loads)
    try_block("Eiffelturm: Referenzergebnis", test_eiffel_reference)
    print("\n========== Projekt 3: Scheibe aus Q4-Elementen ==========")
    try_block("Schritt 1 -- q4/shape.py: Ansatzfunktionen", test_q4_shape_functions)
    try_block("Schritt 2 -- q4/shape.py: Gauss-Punkte", test_q4_gauss)
    try_block("Schritt 3 -- q4/element.py Teil A: Jacobi-Matrix und Gradienten", test_q4_geometry)
    try_block("Schritt 4 -- q4/element.py Teil B: F und E", test_q4_kinematics)
    try_block("Schritt 5 -- q4/material.py", test_material_2d)
    try_block("Schritt 6 -- q4/element.py Teil C: B-Matrix", test_q4_b_matrix)
    try_block("Schritt 7 -- q4/element.py Teil D: f_int und K_e", test_q4_element)
    try_block("Schritt 8 -- q4/assembly.py", test_assembly_q4)
    try_block("Schritt 9 -- q4/mesh.py", test_mesh)
    try_block("Schritt 10 -- newton_raphson.solve und load_stepping.run_general", test_general_solver)
    try_block("Schritt 11 -- Kragbalken: Referenzergebnis", test_cantilever_reference)
    print("\nFertig. Alles OK? Dann committen. Sonst erst den FEHLER beheben.")
