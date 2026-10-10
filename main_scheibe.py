"""
Kragbalken aus Q4-Elementen: links fest eingespannt, rechts nach unten belastet.

Erstes Beispiel mit einem Kontinuumselement statt Staeben. Neu sind Netz
(q4/mesh.py), Ansatzfunktionen (q4/shape.py), Material (q4/material.py), Element
(q4/element.py) und Assemblierung (q4/assembly.py). Newton und Lastschritte
laufen ueber die verallgemeinerten Funktionen newton_raphson.solve und
load_stepping.run_general.
"""
import numpy as np
import matplotlib.pyplot as plt

from q4 import mesh, material, assembly
from loeser import load_stepping
import postprocess

# --- Geometrie & Netz ---
Lx, Ly = 10.0, 1.0                   # Laenge und Hoehe des Balkens
nx, ny = 20, 2                       # Elemente in x- und y-Richtung
t = 1.0                              # Dicke
nodes, elements = mesh.rect_mesh(nx, ny, Lx, Ly)

# --- Material ---
Emod, nu = 1000.0, 0.3
lam, mu = material.lame(Emod, nu)

# --- Randbedingungen: linke Kante (X = 0) fest in x und y ---
n_dofs = 2 * nodes.shape[0]
left_nodes = np.where(np.isclose(nodes[:, 0], 0.0))[0]
fixed_dofs = [d for n in left_nodes for d in (2 * n, 2 * n + 1)]
free_dofs = [d for d in range(n_dofs) if d not in fixed_dofs]

# --- Last: Gesamtkraft P nach unten, gleichmaessig auf die Knoten der rechten Kante ---
P = 1.0                              # mal 0.1 (fast linear) und 3.0 probieren
right_nodes = np.where(np.isclose(nodes[:, 0], Lx))[0]
F_max = np.zeros(n_dofs)
F_max[2 * right_nodes + 1] = -P / len(right_nodes)

# --- Loesen ---
n_steps = 10
assemble_fn = lambda u: assembly.assemble(nodes, elements, u, lam, mu, t)
disp_history, last_residuals = load_stepping.run_general(
    assemble_fn, n_dofs, F_max, free_dofs, n_steps=n_steps
)

# --- Plotten ---
tip = right_nodes[len(right_nodes) // 2]          # Knoten in der Mitte der rechten Kante
load_factor = np.linspace(0.0, 1.0, n_steps + 1)

fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))
axes[0].plot(-disp_history[:, 2 * tip + 1], load_factor * P, "o-")
axes[0].set_xlabel("Durchbiegung am Balkenende (nach unten)")
axes[0].set_ylabel("Last P")
axes[0].set_title("Last-Verschiebungs-Kurve")
axes[0].grid(True)
postprocess.plot_convergence(last_residuals, ax=axes[1])
postprocess.plot_mesh_q4(nodes, elements, disp_history[-1], scale=1.0, ax=axes[2])
plt.tight_layout()
plt.show()

# Zum Vergleich: lineare Balkentheorie (kleine Verformungen, ebener Verzerrungszustand)
w_linear = P * Lx**3 / (3 * (Emod / (1 - nu**2)) * (t * Ly**3 / 12))
print(f"Durchbiegung am Ende: {-disp_history[-1, 2 * tip + 1]:.4f}")
print(f"Lineare Balkentheorie: {w_linear:.4f}")
print(f"Letzter Lastschritt: {len(last_residuals)} Newton-Iterationen.")
