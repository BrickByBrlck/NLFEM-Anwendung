# NLFEM-Anwendung — Von-Mises-Fachwerk

Ein kleines, selbst gebautes Programm, um U1–U3 (Kinematik, Materialgesetz,
Stabelement + Newton-Raphson) an einem konkreten, sichtbaren Beispiel
anzuwenden — auf Wunsch deines Chefs, um zu sehen, "was das eigentlich ist",
nicht als weitere Kursaufgabe. Kein Bezug zu `NLFEM-Lernen`/`NLFEM-Drill`,
kein Blick in `NLFEM-Loesungen` nötig — das hier ist dein eigener Code.

## Was am Ende passiert

Ein Fachwerk aus zwei Stäben, die sich in einer Spitze treffen (klassisches
"Von-Mises-Fachwerk"):

```
        1 (frei, wird belastet)
       / \
      /   \
     0     2   (beide fest eingespannt)
```

Knoten 1 wird schrittweise nach unten belastet. Weil die Stäbe dabei ihren
Winkel ändern, ist die Last-Verschiebungs-Beziehung **nicht linear**, obwohl
das Materialgesetz selbst linear ist (S_tt = E·E_tt) — die Nichtlinearität
steckt komplett in der Geometrie, genau wie in U3. Am Ende siehst du drei
Plots: die Last-Verschiebungs-Kurve, die Konvergenz der letzten
Newton-Iteration, und die verformte Struktur.

## Projektstruktur (so ist ein kleiner FEM-Code normalerweise aufgebaut)

| Datei | Zweck | Musst du ausfüllen? |
|---|---|---|
| `material.py` | Materialgesetz S_tt(E_tt) | Ja (trivial, aber tippen) |
| `element.py` | Stab-Kinematik, Elementresiduum, Elementtangente | **Ja — der Kern** |
| `assembly.py` | lokale → globale Matrizen/Vektoren | Ja |
| `newton_raphson.py` | die Newton-Iteration selbst | **Ja — das, worum es dir ging** |
| `load_stepping.py` | Lastschritte um Newton-Raphson herum | Ja |
| `postprocess.py` | Plots | Nein, fertig — das ist der "Sehen"-Teil |
| `main.py` | Geometrie/Material/Randbedingungen festlegen, alles aufrufen | Nein, fertig |

Diese Aufteilung (Material getrennt von Element, Element getrennt von
Assembly, Assembly getrennt vom Solver) ist kein Selbstzweck — genau diese
Trennung ist auch der Grund, warum man später z.B. nur die Assembly
beschleunigen kann, ohne das Materialgesetz anzufassen.

## Reihenfolge zum Abarbeiten

1. **`material.py`** — zwei Einzeiler (`S_tt = Emod*E_tt` und deren Ableitung).
2. **`element.py`** — die Formeln stehen als Kommentar im Docstring (dieselben
   Größen wie in U3: `L0`, `e`, `d`, `E_tt`, `S_tt`, Elementresiduum,
   `K_mat` + `K_geo`). Übersetze sie in NumPy-Code.
3. **`assembly.py`** — Schleife über Elemente, lokale 4×4/4×1-Beiträge an den
   richtigen globalen Indizes aufaddieren.
4. **`newton_raphson.py`** — die eigentliche Newton-Schleife für **einen**
   Lastschritt: Residuum bilden, Konvergenz prüfen, lineares System lösen,
   Update.
5. **`load_stepping.py`** — Schleife über Lastschritte, die 4. wiederholt
   aufruft und die Historie für die Plots mitschreibt.
6. `python main.py` (bzw. `py main.py`) ausführen.

Solange etwas fehlt, bricht `main.py` mit einer klaren `NotImplementedError`
an genau der Stelle ab, an der du als Nächstes weitermachen musst — das ist
Absicht, kein Bug.

## Selbst-Checks unterwegs (nicht erst am Ende testen)

- **Unverformter Stab:** `element.bar_strain(X1, X2, np.zeros(2), np.zeros(2))`
  muss exakt `0.0` ergeben.
- **Reine Dehnung entlang des Stabs:** Wenn du `u2` parallel zu `e` wählst
  (z.B. `u2 = 0.01 * e`, `u1 = 0`), muss `E_tt` für kleine Werte ungefähr
  `Δ/L0` sein (der quadratische Term ist bei kleiner Dehnung fast egal).
- **Nulllast:** Mit `F_max = 0` muss `load_stepping.run(...)` sofort
  `u = 0` liefern, ohne dass mehr als 1 Newton-Iteration nötig ist.
- **Der eigentliche Korrektheitstest (aus dem Kurs, Punkt 10 in U3):**
  Sobald alles läuft, schau dir den Konvergenz-Plot (mittleres Panel) an.
  `||R||` sollte pro Iteration **quadratisch** fallen — nicht linear. Grob:
  von 1e-2 auf 1e-5 auf 1e-10 in drei Schritten ist quadratisch; braucht es
  zehn Schritte, um von 1e-2 auf 1e-5 zu kommen, ist etwas an der Tangente
  (meistens `K_geo`) falsch. Das ist der zuverlässigste Debug-Trick, den es
  hier gibt — zuverlässiger als jede Formel nochmal nachzurechnen.

## Bonus, wenn es läuft

Erhöhe `F_max_mag` in `main.py` von `15.0` Richtung `35`–`40`. Ab einem
bestimmten Punkt (Snap-Through) wird die Last-Verschiebungs-Kurve nicht mehr
monoton, und die last-gesteuerte Newton-Raphson-Iteration wird bei einem
Lastschritt **nicht mehr konvergieren** — nicht weil dein Code falsch ist,
sondern weil reine Lastkontrolle an einem Grenzpunkt (limit point)
grundsätzlich versagt. Genau dafür gibt es Bogenlängenverfahren
(Arc-Length-Methods), die hier bewusst nicht implementiert sind. Guter
Diskussionspunkt für deinen Chef, falls das Thema noch kommt.

---

# Projekt 2: 2D-Eiffelturm unter Wind

Derselbe Code, eine neue Struktur. `material.py`, `element.py`,
`assembly.py`, `newton_raphson.py` und `load_stepping.py` bleiben
**unveraendert** -- sie funktionieren fuer jedes 2D-Fachwerk, egal ob 3 oder
300 Knoten. Neu sind nur Geometrie und Lasten.

| Datei | Zweck | Musst du ausfuellen? |
|---|---|---|
| `geometry.py` | Knoten, Staebe und Auflager des Turms erzeugen | **Ja** |
| `loads.py` | Windlast und Eigengewicht als Lastvektor | **Ja** |
| `main_eiffel.py` | alles zusammenstecken und plotten | Nein, fertig |

## Reihenfolge

1. `geometry.py`: `half_width` → `eiffel_2d` → `base_dofs`
2. `loads.py`: `wind_load` → `self_weight`
3. Nach jeder Funktion `py selbsttest.py` -- der Abschnitt "Projekt 2" zeigt,
   was schon stimmt
4. `py main_eiffel.py`

Der wichtigste Check ist **"Fachwerk steif?"**: Fehlt ein Stab, kann sich
der Turm ohne Stabdehnung bewegen (Mechanismus), `K_T` wird singulaer und
Newton bricht mit `LinAlgError: Singular matrix` ab. Der Selbsttest sagt dir
dann, wie viele Bewegungsmoeglichkeiten zu viel sind.

## Zum Ausprobieren, wenn es laeuft

- `wind_total` in `main_eiffel.py` auf 20, 50, 100 -- wird die
  Last-Verschiebungs-Kurve krumm?
- `rho_g = 0` (ohne Eigengewicht) vs. `rho_g = 1`
- In `eiffel_2d` testweise eine Diagonale weglassen -> Mechanismus
- `n_levels` auf 3 oder 20

## Moegliche naechste Schritte

- Querschnitt `A` pro Stab (Beine dick, Diagonalen duenn)
- Staebe nach Stabkraft `N` einfaerben (Zug rot, Druck blau)
- Der typische Bogen zwischen den Beinen unten
- 3D-Fachwerk
