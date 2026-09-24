"""Plots -- hier gibt es nichts zu implementieren, das ist der "Sehen"-Teil."""
import numpy as np
import matplotlib.pyplot as plt


def plot_equilibrium_path(load_history, disp_history, dof_index, ax=None):
    """Last-Verschiebungs-Kurve (equilibrium path) an einem DOF."""
    if ax is None:
        _, ax = plt.subplots()
    ax.plot(disp_history[:, dof_index], load_history, "o-")
    ax.set_xlabel("Verschiebung")
    ax.set_ylabel("Last F")
    ax.set_title("Last-Verschiebungs-Kurve")
    ax.grid(True)
    return ax


def plot_convergence(residual_norms, ax=None):
    """Residuumsnorm pro Newton-Iteration, halblogarithmisch.

    Quadratische Konvergenz zeigt sich daran, dass die Kurve nach unten
    zunehmend steiler wird (nicht als Gerade im semilog-Plot).
    """
    if ax is None:
        _, ax = plt.subplots()
    ax.semilogy(range(1, len(residual_norms) + 1), residual_norms, "o-")
    ax.set_xlabel("Newton-Iteration")
    ax.set_ylabel("||R|| (log)")
    ax.set_title("Konvergenz des letzten Lastschritts")
    ax.grid(True, which="both")
    return ax


def plot_structure(nodes, elements, u, scale=1.0, ax=None):
    """Unverformte vs. verformte Struktur."""
    if ax is None:
        _, ax = plt.subplots()
    n_nodes = nodes.shape[0]
    u = np.asarray(u).reshape(n_nodes, 2)
    deformed = nodes + scale * u

    for i, (n1, n2) in enumerate(elements):
        ax.plot(*zip(nodes[n1], nodes[n2]), "k--", alpha=0.5,
                 label="unverformt" if i == 0 else None)
        ax.plot(*zip(deformed[n1], deformed[n2]), "b-o",
                 label="verformt" if i == 0 else None)

    ax.set_aspect("equal")
    ax.legend()
    ax.set_title(f"Struktur (Skalierung x{scale})")
    return ax
