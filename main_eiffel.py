"""
2D-Eiffelturm unter Wind und Eigengewicht.

Verwendet unveraendert dieselben Module wie das Von-Mises-Fachwerk
(element, assembly, newton_raphson, load_stepping) -- neu sind nur die
Geometrie (geometry.py) und die Lasten (loads.py).
"""
import numpy as np
import matplotlib.pyplot as plt

import geometry
import loads
import load_stepping
import postprocess

# --- Geometrie (Einheiten frei gewaehlt, etwa Massstab 1:100) ---
nodes, elements = geometry.eiffel_2d(n_levels=8, height=3.0,
                                     base_half_width=0.625, top_half_width=0.05)

# --- Material & Querschnitt ---
Emod = 1000.0
A = 1.0

# --- Randbedingungen: beide Fusspunkte fest ---
n_dofs = 2 * nodes.shape[0]
fixed_dofs = geometry.base_dofs()
free_dofs = [d for d in range(n_dofs) if d not in fixed_dofs]

# --- Lasten: Wind von links + Eigengewicht ---
wind_total = 10.0                    # gesamte Windkraft; mal 5, 20, 50 probieren
rho_g = 1.0                          # Gewicht pro Volumen (0 = ohne Eigengewicht)
F_max = loads.wind_load(nodes, wind_total) + loads.self_weight(nodes, elements, rho_g, A)

# --- Loesen ---
n_steps = 20
_, disp_history, last_residuals = load_stepping.run(
    nodes, elements, F_max, free_dofs, Emod, A, n_steps=n_steps
)

# --- Plotten ---
apex = nodes.shape[0] - 1            # Spitze = letzter Knoten
load_factor = np.linspace(0.0, 1.0, n_steps + 1)

fig, axes = plt.subplots(1, 3, figsize=(15, 5))
axes[0].plot(disp_history[:, 2 * apex], load_factor * wind_total, "o-")
axes[0].set_xlabel("horizontale Verschiebung der Spitze")
axes[0].set_ylabel("Windkraft gesamt")
axes[0].set_title("Last-Verschiebungs-Kurve (Spitze)")
axes[0].grid(True)
postprocess.plot_convergence(last_residuals, ax=axes[1])
postprocess.plot_structure(nodes, elements, disp_history[-1], scale=1.0, ax=axes[2])
plt.tight_layout()
plt.show()

print(f"Spitze verschiebt sich um {disp_history[-1, 2 * apex]:.4f} (horizontal).")
print(f"Letzter Lastschritt: {len(last_residuals)} Newton-Iterationen.")
