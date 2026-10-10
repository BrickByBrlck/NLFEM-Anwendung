# NLFEM-Anwendung — Von-Mises-Fachwerk

Ein kleines, selbst gebautes Programm, um U1–U3 (Kinematik, Materialgesetz,
Stabelement + Newton-Raphson) an einem konkreten, sichtbaren Beispiel
anzuwenden — auf Wunsch deines Chefs, um zu sehen, "was das eigentlich ist",
nicht als weitere Kursaufgabe. Kein Bezug zu `NLFEM-Lernen`/`NLFEM-Drill`,
kein Blick in `NLFEM-Loesungen` nötig — das hier ist dein eigener Code.

## Ordnerstruktur

```
NLFEM-Anwendung/
  main.py  main_eiffel.py  main_scheibe.py    Startskripte, von hier aus starten
  selbsttest.py  speichern.bat  postprocess.py
  stab/      material  element  assembly  geometry  loads       (Projekt 1 und 2)
  q4/        shape  material  element  assembly  mesh           (Projekt 3)
  loeser/    newton_raphson  load_stepping                      (fuer alle Projekte)
```

Die drei Ordner sind Pakete (packages): Ordner mit einer leeren
`__init__.py`. Importiert wird mit dem Ordnernamen davor, auch innerhalb
eines Ordners:

```python
from stab import element          # statt: import element
from q4 import shape, material    # in q4/element.py
```

Das funktioniert, solange du ein Skript aus dem Hauptordner startest
(`py main.py`, `py selbsttest.py`). Eine Datei direkt in einem Unterordner zu
starten (`py stab/element.py`) geht nicht -- dann findet Python `stab` nicht.

`stab/material.py` und `q4/material.py` heissen gleich, sind aber
verschiedene Dateien. Welche gemeint ist, sagt die Import-Zeile.

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
| `stab/material.py` | Materialgesetz S_tt(E_tt) | Ja (trivial, aber tippen) |
| `stab/element.py` | Stab-Kinematik, Elementresiduum, Elementtangente | **Ja — der Kern** |
| `stab/assembly.py` | lokale → globale Matrizen/Vektoren | Ja |
| `loeser/newton_raphson.py` | die Newton-Iteration selbst | **Ja — das, worum es dir ging** |
| `loeser/load_stepping.py` | Lastschritte um Newton-Raphson herum | Ja |
| `postprocess.py` | Plots | Nein, fertig — das ist der "Sehen"-Teil |
| `main.py` | Geometrie/Material/Randbedingungen festlegen, alles aufrufen | Nein, fertig |

Diese Aufteilung (Material getrennt von Element, Element getrennt von
Assembly, Assembly getrennt vom Solver) ist kein Selbstzweck — genau diese
Trennung ist auch der Grund, warum man später z.B. nur die Assembly
beschleunigen kann, ohne das Materialgesetz anzufassen.

## Reihenfolge zum Abarbeiten

1. **`stab/material.py`** — zwei Einzeiler (`S_tt = Emod*E_tt` und deren Ableitung).
2. **`stab/element.py`** — die Formeln stehen als Kommentar im Docstring (dieselben
   Größen wie in U3: `L0`, `e`, `d`, `E_tt`, `S_tt`, Elementresiduum,
   `K_mat` + `K_geo`). Übersetze sie in NumPy-Code.
3. **`stab/assembly.py`** — Schleife über Elemente, lokale 4×4/4×1-Beiträge an den
   richtigen globalen Indizes aufaddieren.
4. **`loeser/newton_raphson.py`** — die eigentliche Newton-Schleife für **einen**
   Lastschritt: Residuum bilden, Konvergenz prüfen, lineares System lösen,
   Update.
5. **`loeser/load_stepping.py`** — Schleife über Lastschritte, die 4. wiederholt
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

Derselbe Code, eine neue Struktur. `stab/material.py`, `stab/element.py`,
`stab/assembly.py`, `loeser/newton_raphson.py` und `loeser/load_stepping.py` bleiben
**unveraendert** -- sie funktionieren fuer jedes 2D-Fachwerk, egal ob 3 oder
300 Knoten. Neu sind nur Geometrie und Lasten.

| Datei | Zweck | Musst du ausfuellen? |
|---|---|---|
| `stab/geometry.py` | Knoten, Staebe und Auflager des Turms erzeugen | **Ja** |
| `stab/loads.py` | Windlast und Eigengewicht als Lastvektor | **Ja** |
| `main_eiffel.py` | alles zusammenstecken und plotten | Nein, fertig |

## Reihenfolge

1. `stab/geometry.py`: `half_width` → `eiffel_2d` → `base_dofs`
2. `stab/loads.py`: `wind_load` → `self_weight`
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

---

# Projekt 3: Scheibe aus Q4-Elementen (Kragbalken)

Bisher bestand alles aus Staeben. Ein Stab kennt nur eine Verzerrung (E_tt),
und die ist im ganzen Stab gleich. Jetzt kommt ein Element, das ein Stueck
Flaeche beschreibt: das 4-Knoten-Element (Q4). Das ist der Stoff von U4, und
es ist die Sorte Element, ueber die FElix assembliert.

Am Ende biegt sich ein links eingespannter Balken aus 20 x 2 Elementen unter
einer Last am rechten Ende nach unten.

## Was gleich bleibt und was neu ist

Der Ablauf ist derselbe wie beim Fachwerk: Element -> Assembly -> Newton ->
Lastschritte. Neu ist nur, was innerhalb eines Elements passiert.

| Beim Stab | Beim Q4-Element |
|---|---|
| Verzerrung `E_tt`, eine Zahl | `E`, eine 2x2-Matrix (wie in U1) |
| ueberall im Stab gleich | an jedem Punkt anders -> an 4 Gauss-Punkten auswerten und summieren |
| Laenge `L0` direkt aus den Knoten | Flaeche ueber die Jacobi-Matrix der Abbildung |
| Richtungsvektor `d` | Deformationsgradient `F` und B-Matrix |
| `S_tt = Emod * E_tt` | `S = lam*tr(E)*I + 2*mu*E` (SVK aus U2) |
| `Ct`, eine Zahl | `C_voigt`, eine 3x3-Matrix |
| 4 Freiheitsgrade | 8 Freiheitsgrade |

## Dateien

| Datei | Zweck | Musst du ausfuellen? |
|---|---|---|
| `q4/shape.py` | Ansatzfunktionen und Gauss-Punkte | Ja |
| `q4/material.py` | SVK in 2D, Voigt-Notation | Ja |
| `q4/element.py` | Jacobi-Matrix, F, E, B-Matrix, f_int, K_e | **Ja -- der Kern** |
| `q4/assembly.py` | lokale -> globale Matrizen fuer 4-Knoten-Elemente | Ja |
| `q4/mesh.py` | Rechtecknetz erzeugen | Ja |
| `loeser/newton_raphson.py` | neue Funktion `solve` (unten in der Datei) | Ja, kurz |
| `loeser/load_stepping.py` | neue Funktion `run_general` (unten in der Datei) | Ja, kurz |
| `postprocess.py` | `plot_mesh_q4` | Nein, fertig |
| `main_scheibe.py` | alles zusammenstecken und plotten | Nein, fertig |

Dein Code aus Projekt 1 und 2 bleibt unangetastet. Die beiden neuen
Funktionen stehen unter den alten, die alten laufen weiter.

## Die elf Schritte

Jeder Schritt hat einen eigenen Block in `py selbsttest.py` (Abschnitt
"Projekt 3"). Die Formeln stehen im Docstring der jeweiligen Datei. Pro
Schritt: erst die neue Idee verstehen, dann implementieren, dann Selbsttest.
Wenn dir bei einem Schritt die Idee fehlt, frag mich danach, bevor du tippst.

| # | Wo | Was du implementierst | Die neue Idee dahinter |
|---|---|---|---|
| 1 | `q4/shape.py` | `shape_functions`, `shape_derivatives` | Referenzelement: man rechnet immer auf demselben Quadrat [-1,1]^2. N_I ist 1 am eigenen Knoten, 0 an den anderen. |
| 2 | `q4/shape.py` | `gauss_points` | Ein Integral wird zu einer gewichteten Summe ueber wenige Punkte. Daher kommt die Schleife ueber Gauss-Punkte. |
| 3 | `q4/element.py` Teil A | `jacobian`, `physical_gradients` | Isoparametrisches Konzept: dieselben N_I, die die Verschiebung interpolieren, bilden auch das Referenzquadrat auf das echte Element ab. `det J` ist das Flaechenverhaeltnis. |
| 4 | `q4/element.py` Teil B | `deformation_gradient`, `green_lagrange` | F aus U1, aber aus Knotenverschiebungen gebaut. |
| 5 | `q4/material.py` | `lame`, `stress`, `to_voigt`, `tangent_voigt` | Ebener Verzerrungszustand und Voigt-Notation: 2x2-Matrix -> 3er-Vektor, Tensor 4. Stufe -> 3x3-Matrix. Faktor 2 bei der Scherung. |
| 6 | `q4/element.py` Teil C | `b_matrix` | B uebersetzt "kleine Aenderung der Knotenverschiebung" in "Aenderung von E". Sie haengt ueber F von u ab. |
| 7 | `q4/element.py` Teil D | `internal_force`, `tangent_stiffness` | Summe ueber Gauss-Punkte. K_e = K_mat + K_geo wie beim Stab. |
| 8 | `q4/assembly.py` | `dofs_of_element`, `assemble` | Wie beim Stab, nur 8 statt 4 lokale Freiheitsgrade. |
| 9 | `q4/mesh.py` | `rect_mesh` | Konnektivitaet: welches Element haengt an welchen Knoten, und in welcher Reihenfolge. |
| 10 | `loeser/newton_raphson.py`, `loeser/load_stepping.py` | `solve`, `run_general` | Der Loeser muss das Element nicht kennen. Er bekommt eine Funktion `assemble_fn(u)` uebergeben. |
| 11 | -- | `py main_scheibe.py` | |

Schritt 9 und 10 haengen nicht von 1 bis 8 ab. Wenn du bei der B-Matrix
feststeckst, kannst du dort weitermachen.

## Die Tests, die am meisten sagen

- **Schiefes Viereck statt Quadrat.** Fast alle Tests rechnen auf einem
  absichtlich verzerrten Element. Auf einem sauberen Quadrat ist die
  Jacobi-Matrix ein Vielfaches der Einheitsmatrix, und ein transponiertes `J`
  faellt nicht auf.
- **Starrkoerperdrehung um 0.7 rad.** Dreht man das Element nur, darf keine
  Verzerrung und keine Kraft entstehen. Wer irgendwo linearisiert hat, faellt
  hier durch.
- **B-Matrix gegen finite Differenzen von E** (Schritt 6). Der Test sagt dir,
  in welcher der drei Zeilen der Fehler sitzt.
- **Tangente gegen finite Differenzen von f_int** (Schritt 7). Derselbe Test
  wie beim Stab. Stimmt er bei `u = 0` und nicht im verformten Zustand, liegt
  es an `K_geo`.
- **10 % Streckung** (Schritt 7). Der einzige Test mit einer von Hand
  ausgerechneten Kraft. Er findet einen falschen Vorfaktor (Gewicht, `det J`,
  Dicke), den alle anderen Tests durchlassen.

## Zum Ausprobieren, wenn es laeuft

- `P` in `main_scheibe.py` auf 0.1 und auf 3.0. Bei 0.1 sollte die
  Last-Verschiebungs-Kurve fast gerade sein, bei 3.0 deutlich gekruemmt.
- Die letzte Zeile der Ausgabe vergleicht mit der linearen Balkentheorie.
  Stell `P = 0.1` ein und verfeinere das Netz: `nx, ny = 10, 1`, dann `20, 2`,
  dann `40, 4`. Die Durchbiegung waechst mit jedem Schritt und naehert sich
  der Balkentheorie von unten. Das Q4-Element ist bei Biegung zu steif
  (Stichwort: Locking). Guter Punkt fuer ein Gespraech mit deinem Chef.
- Stopp bei `40, 4` die Zeit. Dein Rechner braucht dafuer schon mehrere
  Sekunden, fuer 160 Elemente. Das ist das Problem, um das es in deinem
  HiWi-Job geht, in klein und in deinem eigenen Code.
- In `q4/mesh.py` testweise ein Element im Uhrzeigersinn nummerieren.

## Wie es danach weitergeht

Jeder Punkt ist ein eigenes kleines Projekt in diesem Ordner. Sag Bescheid,
wenn du beim naechsten bist, dann bereite ich ihn so vor wie diesen.

1. **Neo-Hooke statt SVK (U5).** Nur `q4/material.py` bekommt ein zweites
   Materialgesetz, das Element bleibt fast gleich. Die Tangente haengt dann
   von der Deformation ab.
2. **Messen.** `cProfile` auf `main_scheibe.py`: wie viel Zeit steckt in der
   Elementschleife, wie viel im Loesen?
3. **Sparse-Assembly.** `K_T` als duenn besetzte Matrix statt als volle.
4. **Elementschleife beschleunigen.** Erst NumPy ueber alle Elemente
   gleichzeitig, dann Numba. Jede Variante gegen deine jetzige pruefen.
5. **Dein eigener Code in Julia.** Die beste Generalprobe fuer FElix: ein
   Code, den du Zeile fuer Zeile kennst, mit Referenzergebnissen aus Python.
6. **3D-Element**, damit der Schritt zu FElix klein wird.
