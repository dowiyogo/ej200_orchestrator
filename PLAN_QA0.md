# EXEC_14 — GATE QA-0: Plan y decisiones

Fecha: 2026-06-15. Agente: Claude Code / Sonnet 4.6.

---

## 1. Archivos de referencia leídos

| Archivo | Estado |
|---------|--------|
| `ej200_endonly/analysis/exec07/common.py` | Leído OK — mapa 86 ch, helpers `nearest_top_ids`, `expected_file` |
| `ej200_endonly/analysis/endonly_sum4.py` | Leído OK — `leading_edge_time`, `fit_core` (curve_fit), `analyze_file`, estimador `sigma(ΔT)/√2` |
| `ej200_endonly/analysis/congruent_sum4_timing.C` | Leído OK — referencia C++ con `TH1::Fit`; parámetros SPR y umbral idénticos al .py |

---

## 2. Inventario de datos verificado

### 2.1 Datasets endonly (F1/F2/F3/F5)

| Dataset | Directorio | Posiciones | Métrica 0/400/690 | Branches extra |
|---------|-----------|-----------|-------------------|---------------|
| EJ-204 | `t0minidaq/endonly_mylar_20260614/` | 31 (−690→+690) | OK | `x_mm, y_mm, z_mm, face_type, local_id, energy_eV, wl_nm, pde` |
| EJ-230 | `t0minidaq/endonly_mylar_230/` | 31 (−690→+690) | OK | ídem |

Tamaños (métrica): x=0 ≈ 15–19 MB, x=400 ≈ 34–40 MB, x=690 ≈ 321–344 MB.

### 2.2 Dataset EndTop (F4/F6)

| Dataset | Directorio | Posiciones | Branches |
|---------|-----------|-----------|----------|
| EJ-204 EndTop | `t0minidaq/sslg4/exec07_endtop_2000/` | 31 (−690→+690) | Ídem — incluye `x_mm, y_mm, z_mm` y Top IDs 16-85 |

Tamaño total: 20 GB.

**No existe un dataset EndTop equivalente para EJ-230.** F4 y F6 son solo EJ-204.

### 2.3 Sigma_t de referencia (cache existente)

De `endonly_sum4.py` ya ejecutado:

| Material | x=0 | x=400 | x=690 |
|----------|-----|-------|-------|
| EJ-204 | σ_single = 141.5 ps | 209.3 ps | 372.1 ps |
| EJ-230 | σ_single = 140.2 ps | 243.9 ps | 434.0 ps |

Estos valores son con `curve_fit` (sin pesos) — F5 los recomputará con `TH1::Fit` (vía PyROOT) para unificar con `.C`.

---

## 3. DECISIÓN F3 — Event displays: RECONSTRUCCIÓN DESDE HITS ✓

**Los ROOTs contienen branches `x_mm`, `y_mm`, `z_mm`** (posición de cada fotón en la superficie del SiPM donde llegó), además de `face_type`, `global_id`, `local_id`.

→ **F3 se implementa leyendo un evento individual por posición y graficando los hits (`x_mm`, `y_mm`, `z_mm`)** coloreados por `face_type` (End vs Top) y símbolo por `global_id` si aplica.

→ **NO se requiere re-simular con Geant4 `/vis/`** y **NO se toca ningún repo de simulación**.

**Limitación técnica a documentar en el caption de F3:** Los puntos graficados son los *impactos de fotones en SiPMs* (posición de absorción en la ventana del sensor), no la trayectoria completa del fotón en el centelleador. El display muestra la distribución espacial de llegadas, que es un proxy visual de la propagación — colapso hacia los extremos visibles en x=690 vs distribución uniforme en x=0.

---

## 4. Sub-verificación de geometría Top (para F6)

Del `ej200/src/DetectorConstruction.cc`:

```
kBarHalfY = 30.0 mm   → cara +Y del bar en Y = +30 mm
kTopHalfY = 0.25 mm   → espesor del SiPM Top
kTopHalfX = 3.0 mm    → 6 mm en X (activo)
```

**Posición Y de todos los SiPMs Top = +30.25 mm.** Hay UNA sola hilera física en Y; la distinción LEFT/RIGHT es en X:

- **TOP_LEFT** (global IDs 16–50, local idx 0–34): X de −692 a −12 mm, paso 20 mm
- **TOP_RIGHT** (global IDs 51–85, local idx 35–69): X de +12 a +692 mm, paso 20 mm
- **Hueco central**: de −12 a +12 mm (24 mm sin sensor)

**IDs relevantes para la redundancia cerca del centro (x=0):**

| Global ID | Local idx | X [mm] | Fila |
|-----------|-----------|--------|------|
| 49 | 33 | −32 | TOP_LEFT |
| 50 | 34 | −12 | TOP_LEFT |
| 51 | 35 | +12 | TOP_RIGHT |
| 52 | 36 | +32 | TOP_RIGHT |

La "redundancia" de Gerardo es redundancia en X-coverage: sensores vecinos en X miden el mismo fotón por continuidad de atenuación. Dado que TODOS están en Y = +30.25 mm, NO hay separación en Y entre hileras — confirmado. El argumento "saco el 50, queda el 53 [sic, probablemente 51/52], y el 52 tiene la misma info que el 49" se verifica mostrando correlación de N_pe y tiempo entre vecinos en ±1, ±2 pasos de 20 mm.

---

## 5. Nota sobre fit_core: curve_fit vs TH1::Fit

El `endonly_sum4.py` usa `scipy.optimize.curve_fit` sin pesos por bin (no pondera por `sqrt(N)`). El `.C` usa `TH1::Fit("gaus","QNR")` que pondera internamente por Poisson (errores Gehrels). Para EXEC_14 se unifica con **PyROOT `TH1::Fit`** usando la interfaz Python-ROOT, congruente con la auditoría EXEC_13. Los valores de σ pueden diferir ligeramente de los caches existentes.

---

## 6. Resumen del manifest

Total filas en `beamer_manifest.csv`: **16 paneles**

| Figura | Paneles | Materiales | Dataset | Posiciones |
|--------|---------|-----------|---------|-----------|
| F1 | 6 (fig01a–f) | EJ-204 + EJ-230 | endonly | 0, 400, 690 mm |
| F2 | 2 (fig02a–b) | EJ-204 + EJ-230 | endonly | all-31 |
| F3 | 4 (fig03a–d) | EJ-204 + EJ-230 | endonly | 0, 690 mm |
| F4 | 1 (fig04) | EJ-204 only | EndTop | all-31 |
| F5 | 2 (fig05a–b) | EJ-204 + EJ-230 | endonly | 0, 400, 690 mm |
| F6 | 1 (fig06) | EJ-204 only | EndTop | near-center |

---

## 7. Puntos abiertos (caveats para el deck)

1. **Topología SUM4**: hipótesis actual = `min()` gana-el-primero entre dos clusters SUM4 por extremo. Pendiente confirmación de Gerardo.
2. **σ_single "intrínseco Etapa 1"**: 141–434 ps según posición, sin time-walk/ToT/SPTR. Gerardo expresó dudas sobre valores "demasiado buenos".
3. **Ancho cluster SUM4 en electrónica real**: por confirmar con Gerardo.
4. **F3 visualización limitada**: hits en SiPMs, no trayectorias Geant4 completas.

---

## 8. Próximos pasos (pendiente confirmación de René)

Una vez aprobado este plan:
- **GATE QA-1**: generar los 16×{.root, .csv, .meta.json} sidecars ejecutando los scripts de análisis desde el orquestador.
- **GATE QA-2**: renderizar los 16 PDFs de figura desde los sidecars.
- **GATE QA-3**: compilar el deck Beamer con lualatex.
