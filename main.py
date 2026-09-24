"""
Von-Mises-Fachwerk (2 Staebe, ein freier Knoten) -- kleinstmoegliches
Beispiel fuer geometrisch nichtlineares Verhalten mit einem
"Snap-Through"-Punkt.

Knoten 0 und 2 sind Auflager (fest), Knoten 1 ist frei und wird mit einer
nach unten wachsenden Kraft belastet.

Erst main.py ausfuehren, wenn material.py, element.py, assembly.py,
newton_raphson.py und load_stepping.py implementiert sind -- vorher bricht
es mit NotImplementedError an der Stelle ab, wo als Naechstes etwas fehlt.
Das ist Absicht, siehe README.md.
"""
import numpy as np
import matplotlib.pyplot as plt

import load_stepping
import postprocess

# --- Geometrie ---
nodes = np.array([
    [0.0, 0.0],   # Knoten 0: Auflager links
    [1.0, 14],   # Knoten 1: Spitze, frei
    [2.0, 0.0],   # Knoten 2: Auflager rechts
])
elements = [(0, 1), (1, 2),(2,0)]

# --- Material & Querschnitt ---
Emod = 1000.0
A = 1.0

# --- Randbedingungen ---
n_dofs = 2 * nodes.shape[0]
fixed_dofs = [0, 1, 4, 5]           # Knoten 0 (ux,uy) und Knoten 2 (ux,uy) fest
free_dofs = [d for d in range(n_dofs) if d not in fixed_dofs]

# --- Last: nach unten an Knoten 1 (DOF 3 = uy von Knoten 1) ---
F_max_mag = 45                   # sicher unterhalb des Snap-Through-Punkts (~36)
F_max = np.zeros(n_dofs)
F_max[3] = -F_max_mag

# --- Loesen ---
load_history, disp_history, last_residuals = load_stepping.run(
    nodes, elements, F_max, free_dofs, Emod, A, n_steps=20
)

# --- Plotten ---
fig, axes = plt.subplots(1, 3, figsize=(15, 4))
postprocess.plot_equilibrium_path(load_history, disp_history, dof_index=3, ax=axes[0])
postprocess.plot_convergence(last_residuals, ax=axes[1])
postprocess.plot_structure(nodes, elements, disp_history[-1], scale=1.0, ax=axes[2])
plt.tight_layout()
plt.show()

print(f"Letzter Lastschritt: {len(last_residuals)} Newton-Iterationen bis Konvergenz.")
print("Faellt ||R|| im mittleren Plot sichtbar quadratisch? Dann ist K_e richtig.")
