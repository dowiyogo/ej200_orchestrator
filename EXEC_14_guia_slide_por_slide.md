# Guía diapositiva-por-diapositiva — Deck EXEC_14
## "Optical Simulation of Scintillator Bars: EJ-204 vs EJ-230 for the SHiP Timing Detector"

Documento de referencia para René: qué muestra cada diapositiva, qué cantidad física representa, y cómo se calculó. Pensado para que puedas defender cada número ante Gerardo o la colaboración, y para que recuerdes las decisiones de análisis detrás de cada figura. 16 diapositivas.

Convenciones que aplican a todo el deck:
- **Cero temporal**: `time_ns` de cada fotón = `GetGlobalTime()` de Geant4, medido desde la generación del muón primario (verificado en `SiPMSD.cc:80`). Incluye tránsito del muón (~ps), emisión de centelleo, propagación óptica y un jitter por-hit de 20 ps.
- **Estimador de resolución**: σ_t = σ(ΔT_LR)/√2, donde ΔT_LR = t_izquierda − t_derecha entre los dos extremos. El √2 asume dos extremos simétricos e independientes.
- **Todo es Etapa 1 (intrínseca)**: sin time-walk, sin corte ToT, sin SPTR del SiPM. El único término electrónico es el jitter de 20 ps/hit ya dentro de `time_ns`.
- **Materiales**: EJ-204 (azul, OPSC-101, decay datasheet 1.8 ns) y EJ-230 (rojo). Ambos con SSLG4.

---

## Diapositiva 1 — Portada
Título, autores (René Ríos, Gerardo Vásquez), grupo y fecha. Sin contenido técnico.

---

## Diapositiva 2 — Contenido
Índice de las 9 secciones. Cada figura F1–F7 tiene su sección; QA-1c es la verificación de time-of-flight; el resumen cierra.

---

## Diapositiva 3 — Geometría de la barra y readout SiPM
**Qué muestra:** los parámetros de la simulación. No es un resultado, es el setup.

**Contenido y de dónde sale:**
- Barra de 1400 mm (±700 mm en x), recubrimiento Mylar reflectivo. Definido en `DetectorConstruction.cc`.
- v_eff ≈ 277 mm/ns nominal axial. **Importante:** la línea secundaria dice que el onset del ToF mide v_eff ≈ 292 mm/ns (+5%), porque el onset lo marcan los rayos más directos, no el camino promedio. Este 292 sale de la medición de QA-1c (slide 7), no es un supuesto.
- 16 SiPMs End (IDs 0–15: LEFT 0–7, RIGHT 8–15) para timing. 70 SiPMs Top (IDs 16–85, una sola hilera en Y=+30.25 mm, paso 20 mm, gap central de 24 mm) para veto/tracking. Mapa verificado en `exec07/common.py`.
- SUM4 = mínimo sobre clusters de 4 vecinos (topología "gana-el-primero", **pendiente confirmar con Gerardo**).
- Estadística: 3 datasets (EJ-204 End, EJ-230 End, EJ-204 EndTop), 2000 eventos/posición, 31 posiciones (0, ±50…±690 mm).

---

## Diapositiva 4 — F1: perfil de tiempo de llegada, posición central x=0
**Qué muestra:** la distribución del tiempo de llegada de los fotones a los SiPMs End, en el centro de la barra, para ambos materiales. Cuatro paneles: EJ-204 (izq) y EJ-230 (der), cada uno en ventana 0–2 ns (bin 2 ps) y 0–10 ns (bin 10 ps), eje Y logarítmico.

**Cantidad graficada:** t_rel = t_fotón − min(t_End del evento). Se resta el primer fotón de cada evento para quitar el offset de propagación común (~2.8 ns en x=0) y dejar ver la forma de la emisión.

**Cómo se calculó:** por cada evento, se toma el tiempo de cada hit en los End SiPMs, se resta el mínimo del evento, y se acumula en el histograma. El ajuste es un modelo de 1 modo: N·(1−e^(−t/τ_r))·e^(−t/τ_f), ajustado con TH1::Fit (pesos Poisson).

**Resultado y lectura:** EJ-204 ajusta bien (χ²/ndf=1.93, τ_r=1.50 ns, τ_f=1.70 ns). EJ-230 tiene χ²/ndf=17.8 — el ajuste de 1 modo no describe bien su cola. **Punto honesto del caption:** los τ ajustados son del perfil de *llegada* (emisión convolucionada con propagación), NO las constantes de emisión del datasheet (EJ-204 nominal: τ_r=0.7, τ_d=1.8 ns). Por eso τ_r sale 1.5 y no 0.7. Y la cola de EJ-230 en x=0 no puede ser ToF (los extremos son equidistantes en el centro); su origen queda sin determinar — modo lento real o reflexiones. Este es el único candidato vivo a "segundo modo de centelleo" del deck.

---

## Diapositiva 5 — F1: perfil de tiempo de llegada, borde x=+690 mm
**Qué muestra:** lo mismo que el slide 4 pero en el borde de la barra. Aparece una segunda componente clara a ~4.7 ns en ambos materiales.

**Lectura:** el ajuste de 1 modo FALLA aquí (χ²/ndf=886 EJ-204, 721 EJ-230) porque hay estructura que el modelo no captura: el bump a ~4.7 ns. El caption dice explícitamente que **no es centelleo lento** — es time-of-flight, demostrado en el slide siguiente. La línea vertical marca el onset del ToF a ~4.72 ns. La cola entre 3–10 ns es ~98% fotones del extremo cercano dispersados + ~2% del ToF del extremo lejano.

---

## Diapositiva 6 — Diagrama explicativo del ToF
**Qué muestra:** un esquema (no datos) de por qué aparece el bump del slide 5. Diagrama TikZ: la barra, el muón pegando cerca del extremo derecho (x=+690), un camino corto al SiPM cercano (llega ~0 ns) y un camino largo de 1380 mm al SiPM lejano (llega ~4.7 ns después).

**Propósito:** convertir "time-of-flight" de jerga a algo evidente. El mensaje clave: el bump tardío es *la misma luz* llegando al extremo opuesto, no una componente lenta de emisión.

---

## Diapositiva 7 — QA-1c: confirmación del ToF
**Qué muestra:** la prueba numérica de que el bump es ToF. Tres paneles: x=+690 para EJ-204 y EJ-230 separando los hits por extremo (BOTH / END_RIGHT cercano / END_LEFT lejano), más x=0 como control.

**Cómo se calculó la prueba:** se histograma t_rel separando los SiPMs del extremo cercano y los del lejano. Si el bump es ToF, debe aparecer SOLO en el extremo lejano y a un tiempo = distancia/v_eff.

**Resultado:** onset medido del extremo lejano 4.72 ns (EJ-204) y 4.74 ns (EJ-230), contra predicción de 1380 mm / 277 mm·ns⁻¹ ≈ 4.98 ns. Diferencia −5%. El extremo cercano no muestra bump (igual que el control x=0). **ToF confirmado.** El −5% no es error: refleja que el onset lo marcan los rayos más directos (v_eff efectiva ~292 mm/ns), no el camino promedio. Esta es la medición de la que sale el 292 del slide 3.

---

## Diapositiva 8 — F3: distribución espacial de impactos de fotones
**Qué muestra:** dónde llegan los fotones a los SiPMs en un evento individual, para x=0 y x=+690. Cuatro paneles: vista longitudinal (x vs y) y sección transversal (z vs y), para cada posición.

**Cómo se calculó:** se toma un evento representativo y se grafican las coordenadas x_mm/y_mm/z_mm de cada hit (reconstruidas desde el ROOT, sin re-simular).

**Lectura:** en x=0 la distribución es simétrica L/R (480 hits izq, 123 der). En x=+690 colapsa hacia el extremo derecho (45 izq, 8690 der) — visualiza el "colapso por scattering" que pidió Gerardo. **Caveat importante del caption:** los puntos son impactos en la ventana del SiPM, NO trayectorias dentro del volumen. Es un proxy visual de la propagación. (Pendiente confirmar con Gerardo si quería esto o trazas ópticas en el volumen, que necesitarían `/vis/`.)

---

## Diapositiva 9 — F2: σ_t vs N_pe (END SiPMs)
**Qué muestra:** cómo depende la resolución temporal del número de fotones recibidos. Cada punto es UNA posición del muón (31 posiciones). Dos materiales, dos series por material.

**MUY IMPORTANTE — qué NO es:** F2 NO es una suma de fotones sobre el Top (eso es F4). F2 usa los SiPMs **End** (timing). Cada punto: eje X = N_pe que recibió el extremo más débil en esa posición; eje Y = la σ_t de esa posición (el mismo número de F5).

**Dos series:** "weak end" = min(N_pe_izq, N_pe_der), el extremo que limita la resolución; y "mean" = promedio de los dos extremos.

**Cómo se calculó:** σ_t por posición se toma del cálculo de F5 (sigma_t_sum4.csv). El N_pe por extremo se cuenta de los ROOT. Se ajusta el modelo σ_t = √(a²/N_pe + b²), donde a²/N_pe es el término estadístico de conteo de fotones y b sería un piso de dispersión óptica que no mejora con más luz.

**Resultado clave:** a=1312 ps (EJ-204), 1198 ps (EJ-230), **b=0 en ambos**. Es decir: en el rango accesible (N_pe del extremo débil entre 11 y 77), la resolución es **régimen Poisson puro** (σ_t ∝ 1/√N_pe), sin piso óptico visible. **Conclusión:** la resolución en los bordes está limitada por la *poca luz* que llega al extremo lejano, NO por dispersión de propagación. Implicación de diseño: más luz al extremo débil (acoplamiento, más SiPMs, menos atenuación) siempre mejora — no hay saturación en el rango de operación.

**Nota para defender:** la curva de ajuste se extiende hasta N_pe~2000 pero los puntos del extremo débil solo llegan a ~77; más allá es extrapolación del modelo. La serie "mean" sube hacia la derecha porque mezcla posiciones de muy distinto ΔT, no porque la resolución empeore con luz.

**Reconciliación con F7:** la dispersión óptica EXISTE (F7 la mide) pero es subdominante frente a la estadística de fotones en la σ_t del SUM4 (b≈0 aquí).

---

## Diapositiva 10 — F4: perfil de luz del Top vs posición
**Qué muestra:** cuánta luz reciben los SiPMs Top según dónde pega el muón. Es un PERFIL DE LUZ, no de timing (lo dice el título).

**MUY IMPORTANTE — qué es vs F2:** F4 usa los SiPMs **Top** (IDs 16–85). Eje Y = fotoelectrones sumados sobre sensores Top. F2 era timing con los End; F4 es luz con el Top. Son figuras distintas.

**Cómo se calculó:** por posición, se promedia sobre eventos el N_pe total del Top (todos los SiPMs Top) y el N_pe del SiPM más cercano al impacto (`nearest_top_ids`).

**Resultado:** el N_pe total decae suavemente del centro (~4476) al borde (~2658) pero nunca baja de ~2600. Por eso los cortes T4 (≥4 PE) y T20 (≥20 PE) están saturados al 100% en todas partes — el Top siempre tiene mucha más luz que cualquiera de esos umbrales. **Mensaje:** la señal Top es robusta en toda la barra (justifica que el Top no necesita estar tan cerca). La curva "Nearest" oscila por el paso discreto de 20 mm de los sensores. **Pendiente Gerardo:** si T4/T20 era un corte en N_pe (lo implementado) o una eficiencia de trigger SUM4.

---

## Diapositiva 11 — F5: resolución SUM4 L/R, EJ-204
**Qué muestra:** la métrica central del deck. Histogramas de ΔT_LR en x=0, +400, +690 mm, más un gráfico de barras con la σ_t en los tres puntos.

**Cómo se calculó:** ΔT_LR = t_izq − t_der por evento (cada extremo usa el SUM4 = mínimo de sus clusters de 4). Se ajusta el core con gaussiana seeded-at-peak (mediana + MAD robusto, estilo Constanza, vía TH1::Fit). σ_single = σ(ΔT_LR)/√2.

**Resultado:** 141 ps (x=0), 206 ps (x=400), 365 ps (x=690). La resolución se degrada hacia el borde. La eficiencia (fracción de eventos con ambos extremos disparando) cae a 95.5% en x=690. **Etiqueta obligatoria:** intrínseco + jitter 20 ps/hit (Etapa 1), sin time-walk ni readout jitter adicional.

---

## Diapositiva 12 — F5: resolución SUM4 L/R, EJ-230
**Qué muestra:** lo mismo que el slide 11 pero para EJ-230.

**Resultado:** 139 ps (x=0), 247 ps (x=400), 419 ps (x=690). Mismo formato y caveats. La degradación posicional es más fuerte que en EJ-204, y la eficiencia en x=690 cae más (78.5%) — EJ-230 pierde más eventos en el borde.

---

## Diapositiva 13 — Tabla comparativa EJ-204 vs EJ-230
**Qué muestra:** los seis números de σ_t en una tabla, con observaciones.

**Lectura:** resolución central equivalente (141 ≈ 139 ps); EJ-204 ~15% mejor en el borde (365 vs 419 ps); EJ-230 degrada más rápido con la posición (+201% vs +159% del centro al borde). Recordatorio de que los τ_f son del perfil de llegada, no del datasheet. Caveats de Etapa 1 listados: posible doble conteo del jitter, sin time-walk, sin readout jitter, topología SUM4 pendiente.

**Para defender:** estos números NO son la resolución del detector real. Falta SPTR (~106 ps), electrónica real, time-walk. Son límites intrínsecos de simulación. Preséntalos siempre con esa etiqueta.

---

## Diapositiva 14 — F6: redundancia de SiPMs Top
**Qué muestra:** correlación de N_pe entre pares de SiPMs Top vecinos, para ver cuáles dan información redundante. Tres scatter 2D (mapas de densidad).

**Cómo se calculó:** por evento, se cuenta N_pe de cada SiPM del par, y se calcula el coeficiente de correlación de Pearson entre los dos.

**Resultado y lectura — esto corrigió la intuición inicial:**
- Par A (IDs 49/52, simétricos al centro, lados opuestos): r=−0.37, **anticorrelados/complementarios** — compiten por los mismos fotones (un fotón al 49 es uno que no llegó al 52).
- Par B (IDs 50/51, cruzando el gap central): r=+0.88, redundantes.
- Par C (IDs 47/49, mismo lado): r=+0.90, redundantes (control).

**Mensaje:** la redundancia real está entre vecinos del MISMO lado (se puede prescindir de uno), no entre los simétricos al centro. Esto matiza lo que Gerardo dijo de memoria (que 49 y 52 daban la misma info): los datos dicen que son complementarios, no redundantes.

---

## Diapositiva 15 — F7: tiempo de llegada por order-statistic, x=−690
**Qué muestra:** para los 9 SiPMs Top más cercanos al impacto (rejilla 3×3), el tiempo medio de llegada del n-ésimo fotón en función de n, comparado con un piso estadístico teórico.

**Cómo se calculó:** por canal, se ordenan los tiempos de llegada de cada evento y se promedia el n-ésimo (⟨t_n⟩). Las barras de error son la RMS evento-a-evento. La banda punteada es el piso estadístico σ_stat(t_n) = √n·τ_d/⟨N_pe⟩ con τ_d=1.8 ns (el límite si los tiempos fueran puras order statistics de la emisión, sin dispersión óptica).

**Lectura — el corazón físico:** la brecha entre la RMS observada y el piso estadístico AÍSLA la contribución de dispersión de camino óptico. El piso es el límite ideal; los puntos incluyen propagación; la brecha es la física de propagación. Es la misma dispersión que ensancha σ_t y que mueve v_eff del onset respecto al nominal. En x=−690 el cluster usa IDs 16–24 (truncado en el borde izquierdo, no hay canales más allá).

**Reconciliación con F2:** la dispersión EXISTE aquí (RMS > piso) pero es subdominante a la estadística Poisson en la σ_t del SUM4 (b≈0 en F2). F7 muestra que el efecto está presente; F2 muestra que aún no domina.

---

## Diapositiva 16 — Resumen y preguntas abiertas
**Resultados principales (5):**
1. EJ-204 ≈ EJ-230 en el centro; EJ-204 mejor en bordes.
2. La componente de F1 a ~4.7 ns en el borde es ToF geométrico, no centelleo lento.
3. v_eff nominal 277 mm/ns; el onset da 292 mm/ns (rayos más directos).
4. F7: la RMS de ⟨t_n⟩ excede el piso estadístico — cuantifica la dispersión óptica, misma física que la degradación de σ_t en bordes.
5. Top robusto (>20 PE); sensores adyacentes al gap muy correlacionados.

**Pendientes — confirmación de Gerardo (6):** topología SUM4; definición T4/T20; jitter 20 ps simple o doble; F3 impactos vs trazas; anticorrelación del Par A (geometría o reflexiones); cola de EJ-230 en x=0 (modo lento o reflexiones).

**Próximos pasos:** dataset EJ-230 EndTop; Etapa 2 (time-walk, readout jitter, electrónica real); deconvolución para recuperar constantes de emisión.

---

## Apéndice — hilo argumental del deck (para tu presentación)

El deck cuenta una historia con un hilo central que conviene resaltar al presentar: **la resolución temporal en los bordes de la barra está limitada por la cantidad de luz que llega al extremo lejano, no por la dispersión de propagación.** Cómo lo construye cada figura:

- F5 establece el hecho: σ_t se degrada del centro al borde (141→365 ps EJ-204).
- F1 + QA-1c descartan una explicación falsa: la estructura tardía en el borde es geometría (ToF), no física de centelleo.
- F7 mide que la dispersión óptica existe (brecha sobre el piso estadístico).
- F2 cierra el argumento: esa dispersión es subdominante (b≈0); lo que limita es la estadística de fotones (Poisson). Más luz siempre ayuda.
- F3 lo visualiza: en el borde, casi toda la luz va a un solo extremo; el otro queda con poca.
- F4 y F6 son sobre el Top (diseño del detector), complementan pero no son parte del hilo de timing.

Si presentas en ese orden lógico, cada figura responde una pregunta que abre la anterior, y la conclusión (más luz al extremo débil mejora el timing en bordes) sale natural y es accionable para el diseño.

Todos los números del deck son intrínsecos de Etapa 1. La resolución del detector real será peor (SPTR ~106 ps domina). Preséntalo siempre con esa salvedad.
