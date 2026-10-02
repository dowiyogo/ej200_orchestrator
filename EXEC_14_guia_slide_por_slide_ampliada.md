# Guía ampliada, diapositiva por diapositiva — Deck EXEC_14
## *Optical Simulation of Scintillator Bars: EJ-204 vs EJ-230 for the SHiP Timing Detector*

Esta versión está escrita para poder leer y defender las figuras sin depender de abreviaturas ambiguas. En particular, aclara exactamente qué significa **END** en F1 y F2, qué extremo es cercano o lejano en cada posición, qué representa cada punto y qué conclusiones sí —y cuáles no— se pueden extraer.

---

# 1. Respuesta directa a la duda sobre “END”

## Diapositiva 4: ¿qué END se usa?

En la diapositiva 4, **END significa los dos extremos físicos juntos**:

- `END_LEFT`: SiPM IDs 0–7, en `x = −700 mm`.
- `END_RIGHT`: SiPM IDs 8–15, en `x = +700 mm`.
- `all End SiPMs`: unión de ambos conjuntos, IDs 0–15.

Por lo tanto, en F1 **no se muestra END_LEFT por separado de END_RIGHT**. Para cada material, se mezclan en un único histograma todos los fotones detectados en cualquiera de los dos extremos.

La frase del título

> EJ-204 (left) vs EJ-230 (right)

no se refiere a los extremos físicos. Significa solamente:

- EJ-204 está en la **columna izquierda de la diapositiva**.
- EJ-230 está en la **columna derecha de la diapositiva**.

La definición correcta del cero temporal es

\[
t_{0,\mathrm{END}}^{(e)}
=
\min_{h\in e,\;\mathrm{ID}(h)\in[0,15]} t_h,
\]

es decir, para cada evento `e` se busca **el primer fotón detectado entre ambos extremos**. Luego cada hit de END_LEFT y END_RIGHT entra al histograma mediante

\[
t_{\mathrm{rel},h}=t_h-t_{0,\mathrm{END}}^{(e)}.
\]

En `x = 0` los extremos son geométricamente equivalentes. El primer fotón puede pertenecer a END_LEFT o END_RIGHT según la fluctuación del evento.

## Diapositiva 9: ¿qué END se usa?

En la diapositiva 9 también se usan **ambos extremos**, pero aquí no se mezclan inmediatamente. Primero se calculan por separado:

\[
\langle N_{pe}^{L}\rangle_x
\quad\text{para END_LEFT, IDs 0–7,}
\]

\[
\langle N_{pe}^{R}\rangle_x
\quad\text{para END_RIGHT, IDs 8–15.}
\]

Después se construyen dos variables horizontales distintas:

\[
N_{pe}^{\mathrm{weak}}
=
\min\!\left(
\langle N_{pe}^{L}\rangle_x,
\langle N_{pe}^{R}\rangle_x
\right),
\]

\[
N_{pe}^{\mathrm{mean}}
=
\frac{
\langle N_{pe}^{L}\rangle_x+
\langle N_{pe}^{R}\rangle_x
}{2}.
\]

Cada posición del muón produce la misma resolución temporal `σ_t`, pero aparece **dos veces** en la figura:

- marcador oscuro: `σ_t` contra `Npe^weak`; esta es la serie físicamente relevante para el ajuste;
- marcador tenue: la misma `σ_t` contra `Npe^mean`; es una referencia para mostrar que el promedio de luz puede ser engañoso.

El extremo débil coincide aproximadamente con el extremo lejano:

| Posición del muón | Extremo cercano | Extremo lejano / débil |
|---|---|---|
| `x > 0` | END_RIGHT | END_LEFT |
| `x < 0` | END_LEFT | END_RIGHT |
| `x = 0` | ninguno; son equivalentes | ninguno; son equivalentes |

---

# 2. Convenciones globales del deck

## 2.1 Coordenada longitudinal

La barra se extiende a lo largo de `x`:

- extremo izquierdo: `x = −700 mm`;
- centro: `x = 0`;
- extremo derecho: `x = +700 mm`.

Una inyección en `x = +690 mm` ocurre a unos 10 mm del extremo derecho y a aproximadamente 1380 mm del extremo izquierdo.

## 2.2 Grupos de SiPM

| Nombre | IDs | Ubicación | Función principal |
|---|---:|---|---|
| END_LEFT | 0–7 | cara izquierda de la barra | timing |
| END_RIGHT | 8–15 | cara derecha de la barra | timing |
| END / all End SiPMs | 0–15 | ambos extremos | término colectivo |
| TOP_LEFT | 16–50 | fila superior, lado izquierdo | veto/tracking/luz |
| TOP_RIGHT | 51–85 | fila superior, lado derecho | veto/tracking/luz |
| TOP | 16–85 | todos los sensores superiores | término colectivo |

**Regla de lectura:** cuando una gráfica dice solamente `END`, no significa automáticamente izquierda o derecha. Debe leerse como el conjunto definido en la leyenda de esa figura.

## 2.3 Tiempo almacenado

El `time_ns` de cada hit corresponde al tiempo global de Geant4 desde la generación del primario. Incluye, de acuerdo con la configuración del análisis:

- propagación del muón hasta la deposición;
- tiempo de emisión del fotón de centelleo;
- propagación óptica;
- reflexiones y trayectorias indirectas;
- jitter por hit de 20 ps incluido en esta etapa.

No incluye todavía una cadena electrónica real completa, corrección de time-walk, SPTR real del SiPM ni jitter adicional de lectura.

## 2.4 Resolución temporal usada

Se define

\[
\Delta T_{LR}=T_L-T_R,
\]

con `T_L` construido desde END_LEFT y `T_R` desde END_RIGHT. Luego

\[
\sigma_t=\frac{\sigma(\Delta T_{LR})}{\sqrt{2}}.
\]

El factor `1/√2` es exacto solamente si ambos extremos tienen resoluciones iguales e independientes. Cerca del borde esa simetría deja de ser buena, porque el extremo lejano recibe mucha menos luz. Por eso, en los bordes, `σ_t` debe interpretarse como una **métrica equivalente o convencional de un extremo**, no como la resolución real individual de END_LEFT y END_RIGHT por separado.

## 2.5 Qué significa `Npe`

En estas simulaciones `Npe` es el número de fotones/hits detectados usado como proxy de fotoelectrones. Según la figura puede significar:

- suma en END_LEFT;
- suma en END_RIGHT;
- mínimo entre ambos extremos;
- promedio entre ambos extremos;
- suma de todos los TOP;
- cuenta de un único SiPM TOP.

Siempre hay que leer la leyenda antes de interpretar el eje.

---

# 3. Guía diapositiva por diapositiva

## Diapositiva 1 — Portada

No contiene resultados. Presenta el propósito general: comparar la respuesta óptica y temporal simulada de barras EJ-204 y EJ-230 para el detector de timing de SHiP.

**Frase para presentarla:** “Esta comparación se concentra en perfiles de llegada, propagación geométrica, producción de luz y resolución temporal intrínseca antes de incorporar la electrónica real.”

---

## Diapositiva 2 — Contenido

Ordena el análisis en tres familias:

1. configuración y geometría;
2. propagación y distribución de la luz, F1/F3/F4/F6/F7;
3. rendimiento temporal, F2/F5.

QA-1c no es una figura independiente de rendimiento: es una comprobación específica para interpretar correctamente la segunda componente de F1.

---

## Diapositiva 3 — Geometría de la barra y readout

### Pregunta que responde

¿Qué se simuló y dónde están los sensores?

### Geometría

- longitud activa: 1400 mm;
- eje longitudinal: `x`;
- centro: `x = 0`;
- caras: `x = ±700 mm`;
- recubrimiento reflectivo de Mylar.

### END

Hay 16 sensores de extremo:

- 8 en la cara izquierda;
- 8 en la cara derecha.

Estos son los sensores usados en F1, F2 y F5.

### TOP

Hay 70 sensores en la fila superior. Estos aparecen en F4, F6 y F7. No deben confundirse con los sensores END.

### Velocidad efectiva

El valor `277 mm/ns` es una velocidad efectiva axial nominal o promedio. No es necesariamente la velocidad microscópica de la luz en el material. Resume caminos ópticos que incluyen ángulos y reflexiones.

El valor `292 mm/ns` proviene del onset medido en QA-1c. El onset está fijado por los fotones con trayectorias más directas, por lo que puede corresponder a una velocidad efectiva mayor que la del camino promedio.

### SUM4

El deck lo describe como el mínimo sobre clusters de cuatro vecinos. La topología exacta —si se suman señales, si gana el primer cluster o si existe otro criterio— sigue marcada como pendiente de confirmación.

### Qué no concluir

No se debe interpretar `v_eff = 292 mm/ns` como una nueva constante del material. Es una velocidad efectiva asociada al primer frente de llegada en esta geometría.

---

## Diapositiva 4 — F1: perfil relativo de llegada en `x = 0`

### Pregunta que responde

¿Cómo se distribuyen temporalmente los fotones que llegan a los extremos cuando el muón cruza el centro de la barra?

### Qué significa cada columna

- columna izquierda de la diapositiva: EJ-204;
- columna derecha de la diapositiva: EJ-230.

Esto no significa END_LEFT y END_RIGHT.

### Qué sensores entran

En cada material entran todos los hits de:

\[
\mathrm{END}=\mathrm{END\_LEFT}\cup\mathrm{END\_RIGHT}.
\]

La propia leyenda de los plots dice `data (all End SiPMs)`.

### Definición de `t_rel`

Para cada evento se encuentra el primer hit entre IDs 0–15 y se fija como cero. Después se llena el histograma con todos los hits de ambos extremos.

Consecuencias importantes:

1. el cero no es el tiempo de creación del fotón;
2. el cero no es el tiempo absoluto de llegada;
3. cada evento aporta al menos un hit exactamente en `t_rel = 0`;
4. la resta elimina el desplazamiento común de propagación, pero también introduce un sesgo de estadística de orden porque se referencia al mínimo del evento.

### Paneles superiores

- ventana `0–2 ns`;
- binning de 2 ps;
- eje vertical logarítmico.

Sirven para ampliar el frente de subida y la zona de máximo. El bin muy fino permite observar detalles del onset, pero aumenta las fluctuaciones bin a bin.

### Paneles inferiores

- ventana `0–10 ns`;
- binning de 10 ps;
- eje vertical logarítmico.

Sirven para estudiar la cola tardía. No se debe comparar directamente la altura “counts/bin” entre paneles superiores e inferiores porque el ancho de bin es distinto.

### Curvas

- datos azules: EJ-204;
- datos rojos: EJ-230;
- curva verde: ajuste de un solo modo

\[
f(t)=N\left(1-e^{-t/\tau_r}\right)e^{-t/\tau_f}.
\]

`τ_r` controla la subida aparente y `τ_f` la caída aparente del perfil de llegada.

### Por qué esos `τ` no son los del datasheet

El ajuste no observa solamente la emisión intrínseca. Observa

\[
\text{emisión}\otimes\text{propagación óptica}\otimes\text{selección del primer hit}.
\]

Por eso el `τ_r` ajustado puede ser aproximadamente 1.5 ns aunque el datasheet de EJ-204 indique 0.7 ns.

### Calidad del ajuste

- EJ-204, ventana larga: `χ²/ndf ≈ 1.93`; el modelo describe razonablemente la forma global, aunque no es perfecto.
- EJ-230, ventana larga: `χ²/ndf ≈ 17.8`; el modelo de un solo modo no describe la cola.

Los parámetros también cambian al pasar de la ventana corta a la larga. Eso es otra señal de que el modelo es fenomenológico y no una extracción directa de constantes fundamentales.

### Interpretación física

En `x = 0`, END_LEFT y END_RIGHT están a la misma distancia. Por eso una cola adicional de EJ-230 no puede explicarse como la diferencia de tiempo entre un extremo cercano y otro lejano. Las posibilidades abiertas son, por ejemplo:

- trayectorias ópticas reflejadas;
- una componente lenta efectiva;
- limitaciones del modelo de un solo modo;
- efecto de la referencia al primer fotón.

### Frase para presentarla

“Cada material combina ambos extremos; el cero es el primer fotón detectado en cualquiera de ellos. EJ-204 se aproxima a un perfil de un modo, mientras EJ-230 conserva una cola que el modelo simple no explica.”

---

## Diapositiva 5 — F1 en `x = +690 mm`

### Pregunta que responde

¿Qué cambia cuando el muón cruza muy cerca del extremo derecho?

### Mapeo geométrico

En esta posición:

- END_RIGHT es el extremo cercano;
- END_LEFT es el extremo lejano.

Sin embargo, igual que en la diapositiva 4, F1 mezcla en un solo histograma los fotones de ambos extremos.

### Estructura temporal

Los fotones del extremo cercano pueden llegar casi inmediatamente respecto del primer hit del evento. Los del extremo lejano necesitan atravesar aproximadamente 1380 mm y empiezan a aparecer unos 4.7 ns más tarde.

La segunda componente no necesariamente se observa como un pico separado muy estrecho. En escala logarítmica se manifiesta como una nueva población o cambio de pendiente que comienza cerca del onset marcado.

### Por qué falla el ajuste

El modelo de un modo supone una sola población temporal. En el borde hay al menos dos poblaciones geométricas superpuestas:

1. fotones del extremo cercano, incluidos caminos directos y reflejados;
2. fotones del extremo lejano, retrasados por el recorrido longitudinal.

Los valores enormes de `χ²/ndf` indican que el ajuste no debe usarse para interpretar `τ_r` o `τ_f` como parámetros físicos confiables en esta posición.

### Qué significa `t ≈ 0`

Es cero relativo. El fotón cercano no llega literalmente en tiempo global cero; llega cerca del primer hit del evento después de restar `t0,END`.

### Frase para presentarla

“En el borde, ambos extremos siguen superpuestos en F1; el extremo derecho domina el frente temprano y el izquierdo introduce una población tardía cuyo onset está en aproximadamente 4.7 ns.”

---

## Diapositiva 6 — Esquema geométrico del ToF

### Qué representa

Es un diagrama explicativo, no un histograma de datos.

- el muón cruza en `x = +690 mm`;
- al END_RIGHT le quedan aproximadamente 10 mm;
- al END_LEFT le quedan aproximadamente 1380 mm.

La diferencia esperada es

\[
\Delta t\approx\frac{1380\ \mathrm{mm}}{v_{\mathrm{eff}}}.
\]

Con `v_eff = 277 mm/ns`, se obtiene aproximadamente 4.98 ns.

### Mensaje central

Los fotones tardíos no son necesariamente una segunda emisión del material. Pueden ser los mismos fotones de centelleo detectados en la cara opuesta después de un recorrido mucho mayor.

### Precisión conceptual

La flecha larga representa una distancia longitudinal efectiva simplificada. Los fotones reales siguen trayectorias tridimensionales y pueden reflejarse, por lo que la distribución tiene anchura y el onset puede diferir del valor promedio.

---

## Diapositiva 7 — QA-1c: separación por extremo

### Pregunta que responde

¿La población de 4.7 ns proviene realmente del extremo lejano?

### Significado de las curvas

En los paneles de `x = +690 mm`:

- gris `BOTH`: unión o suma bin a bin de los hits de ambos extremos;
- rojo discontinuo `END_RIGHT`: extremo cercano;
- azul punteado `END_LEFT`: extremo lejano.

`BOTH` no significa una coincidencia temporal ni `T_L − T_R`; significa simplemente que ambas poblaciones de hits se muestran juntas.

### Evidencia

La curva azul del extremo lejano comienza alrededor de 4.72–4.74 ns. La curva roja del extremo cercano ya está presente desde el comienzo.

Eso prueba que la nueva población de F1 coincide con la llegada de la luz al extremo opuesto.

### Panel de control `x = 0`

En el centro no existe un extremo cercano y otro lejano. Las curvas izquierda y derecha deben superponerse aproximadamente. Si la leyenda conserva las palabras `near` y `far`, son etiquetas heredadas de la convención de `x = +690 mm`, no una descripción física del panel central.

### Predicción y medición

\[
\Delta t_{\mathrm{pred}}\approx4.98\ \mathrm{ns},
\]

mientras el onset medido es aproximadamente:

- EJ-204: 4.72 ns;
- EJ-230: 4.74 ns.

### Interpretación del 5 %

El primer fotón del extremo lejano selecciona las trayectorias más rápidas. El valor `292 mm/ns` se obtiene aproximadamente de

\[
v_{\mathrm{eff}}^{\mathrm{onset}}
\approx
\frac{1380\ \mathrm{mm}}{4.72\ \mathrm{ns}}
\approx292\ \mathrm{mm/ns}.
\]

No es contradictorio con el promedio de 277 mm/ns: miden partes distintas de la distribución de caminos.

### Qué no concluir

La coincidencia del onset con ToF no implica que toda la cola tardía sea exclusivamente luz del extremo lejano. El extremo cercano también puede producir fotones tardíos por reflexiones. QA-1c identifica el origen de la **nueva población a partir del onset**, no necesariamente de cada fotón individual de la cola.

---

## Diapositiva 8 — F3: impactos espaciales en los SiPM

### Pregunta que responde

¿En qué ventanas de SiPM terminan los fotones de un evento representativo?

### Qué son los puntos

Cada punto es una coordenada de impacto registrada en la ventana sensible de un SiPM. No es un punto de una trayectoria óptica dentro de la barra.

### Vista longitudinal `x_hit` vs `y_hit`

Los impactos aparecen cerca de las caras `x ≈ −700 mm` y `x ≈ +700 mm`, porque solamente se muestran detecciones en los SiPM de extremo.

- azul: END_LEFT;
- rojo: END_RIGHT;
- línea vertical negra: posición del muón `x_gun`.

### Sección transversal `z_hit` vs `y_hit`

Muestra cómo se distribuyen los impactos dentro del área transversal de las ventanas de los sensores.

### Evento central

El evento escogido en `x = 0` tiene 480 hits en END_LEFT y 123 en END_RIGHT. Por lo tanto, no es correcto describirlo como numéricamente simétrico. Lo correcto es decir:

- el muón está geométricamente centrado;
- ambos extremos reciben luz;
- un único evento puede mostrar una fluctuación considerable entre izquierda y derecha.

Una afirmación de simetría cuantitativa requeriría promediar muchos eventos, no seleccionar uno.

### Evento de borde

En `x = +690 mm` aparecen 8690 hits en END_RIGHT y solamente 45 en END_LEFT. Esta diferencia visualiza por qué el extremo lejano limita el timing.

### Qué no concluir

No se puede reconstruir el camino de reflexión de cada fotón a partir de estos puntos finales. Tampoco se deben usar estos dos eventos para calcular eficiencias o cocientes promedio.

### Frase para presentarla

“F3 no dibuja rayos dentro de la barra; dibuja dónde terminan. En el borde, casi todos los impactos se concentran en END_RIGHT y END_LEFT queda estadísticamente pobre.”

---

## Diapositiva 9 — F2: `σ_t` contra luz en END

### Pregunta que responde

¿La resolución temporal empeora porque existe poca luz o porque hay un piso irreducible de propagación?

### Qué representa cada punto

Cada posición `x_gun` produce una distribución de `ΔT_LR` sobre aproximadamente 2000 eventos. De esa distribución se obtiene una sola `σ_t`.

Para la misma posición también se calcula la luz media en cada extremo. Por eso cada posición genera:

- un valor de `σ_t`;
- un `Npe^weak`;
- un `Npe^mean`.

### Marcadores oscuros: extremo débil

\[
N_{pe}^{\mathrm{weak}}
=
\min(\langle N_{pe}^{L}\rangle,
     \langle N_{pe}^{R}\rangle).
\]

Esta es la variable que sigue la rama descendente y que se usa para la interpretación física principal.

En el borde:

- el extremo cercano puede recibir miles de fotones;
- el extremo lejano puede recibir solamente unas decenas;
- `σ_t` queda controlada por el extremo lejano.

### Marcadores tenues: promedio de ambos extremos

\[
N_{pe}^{\mathrm{mean}}
=
\frac{\langle N_{pe}^{L}\rangle+
\langle N_{pe}^{R}\rangle}{2}.
\]

Son los mismos valores de `σ_t` dibujados en otra coordenada horizontal. La rama tenue aumenta hacia la derecha porque, cerca de un borde, el enorme rendimiento del extremo cercano hace crecer el promedio aunque el extremo lejano siga siendo pobre.

Por eso el promedio no es una buena variable para predecir la resolución L/R.

### Ejemplo conceptual

Supóngase:

\[
N_L=20,\qquad N_R=1000.
\]

Entonces

\[
N_{\mathrm{weak}}=20,
\qquad
N_{\mathrm{mean}}=510.
\]

El promedio sugeriría “mucha luz”, pero el tiempo L/R sigue limitado por los 20 fotones del lado débil.

### Ajuste

El modelo es

\[
\sigma_t(N_{pe})
=
\sqrt{\frac{a^2}{N_{pe}}+b^2}.
\]

- `a/√Npe`: contribución que mejora con estadística de fotones;
- `b`: piso constante que no mejora agregando luz.

La curva verde sigue la serie de `Npe^weak`, no la serie tenue de `Npe^mean`.

### Resultado

El ajuste lleva `b` al límite cero en ambos materiales. La lectura prudente es:

> dentro del rango medido, no se resuelve un piso constante distinto de cero y domina una tendencia aproximadamente proporcional a `1/√Npe`.

No conviene afirmar “Poisson puro” como una verdad exacta, porque:

- `χ²/ndf ≈ 5.2` es alto;
- el modelo simple no describe todos los detalles;
- `b = 0` puede ser el resultado de un parámetro en el borde permitido;
- la curva se dibuja mucho más allá del rango realmente poblado por `Npe^weak`.

### Rango observado y extrapolación

Los puntos del extremo débil cubren aproximadamente `Npe = 11–77`. La parte de la curva verde que se extiende hacia cientos o miles de fotoelectrones es extrapolación, no dato.

### Conclusión de diseño

En el rango simulado, mejorar la luz que llega al extremo lejano debería mejorar la resolución. Esto puede lograrse, por ejemplo, mediante mejor acoplamiento, mayor cobertura o menor pérdida óptica. La figura no demuestra que esta mejora continúe indefinidamente fuera del rango estudiado.

### Relación con F7

F7 muestra que existe dispersión de caminos ópticos. F2 indica que esa dispersión no aparece como un piso `b` claramente resuelto en esta métrica y en este rango. La relación es cualitativa, además de que F7 usa sensores TOP y F2 usa END.

### Frase para presentarla

“Ambos extremos se calculan por separado; la curva relevante usa el que recibe menos luz. El promedio puede ser enorme por el extremo cercano y aun así el timing ser malo.”

---

## Diapositiva 10 — F4: luz en los SiPM TOP

### Pregunta que responde

¿Cuánta luz registra la fila superior según la posición del muón?

### Diferencia con F2

- F2: sensores END, eje vertical temporal `σ_t`;
- F4: sensores TOP, eje vertical de luz `Npe`.

F4 no es una gráfica de resolución temporal.

### Curva total

`Total (all Top SiPMs)` suma la luz de IDs 16–85 y promedia por posición. Se mantiene del orden de varios miles de fotoelectrones.

### Curva nearest

`Nearest SiPM` usa el sensor superior más cercano a la posición del muón. La forma de diente de sierra aparece porque los sensores están separados por un pitch discreto de 20 mm:

- cuando el muón pasa justo bajo un sensor, ese sensor recibe más luz;
- entre dos sensores, el sensor “más cercano” queda más lejos y la cuenta baja.

### Escala logarítmica

El eje vertical es logarítmico. La separación visual entre la curva total y la nearest representa factores multiplicativos, no diferencias lineales.

### T4 y T20

En la implementación mostrada, T4 y T20 son perfiles calculados después de seleccionar eventos con `Npe_total ≥ 4` o `≥ 20`. Como todos los eventos superan ambos umbrales, las curvas quedan exactamente encima de la curva total.

Eso es lo que significa “saturadas al 100 %”: la eficiencia del corte es 100 %, no que el SiPM electrónico esté saturado.

### Limitación

Solo se muestra EJ-204 EndTop. No permite comparar directamente la respuesta TOP de EJ-204 con EJ-230.

---

## Diapositiva 11 — F5: resolución SUM4 de EJ-204

### Pregunta que responde

¿Cuánto varía el tiempo izquierda-derecha evento a evento y cómo cambia esa anchura con la posición?

### Histogramas

El eje horizontal es

\[
\Delta T_{LR}=T_L-T_R.
\]

El eje vertical es el número de eventos.

### Significado del centro del pico

- en `x = 0`, las distancias son similares y el pico queda cerca de cero;
- en `x > 0`, END_RIGHT es cercano y `T_R` llega antes;
- por lo tanto `T_L-T_R` es positivo y el pico se desplaza a valores positivos.

El desplazamiento del pico mide principalmente geometría/posición. No es por sí solo una pérdida de resolución.

### Significado de la anchura

La anchura evento a evento es la que se usa para obtener `σ_t`. Un pico más ancho significa peor precisión temporal.

### Valores

- `x = 0`: 141 ± 4 ps;
- `x = +400 mm`: 206 ± 6 ps;
- `x = +690 mm`: 365 ± 10 ps.

### Colas y no gaussianidad

En el borde aparece una cola hacia tiempos más grandes. El número reportado proviene de un ajuste al core. Por eso resume la parte central de la distribución y no describe completamente los outliers tardíos.

### Eficiencia

La eficiencia indicada es la fracción de eventos en que se pudo construir el estimador con información válida en ambos extremos.

- 100 % en centro y `+400 mm`;
- 95.5 % en `+690 mm`.

La caída ocurre porque el extremo lejano puede no reunir suficientes hits para el criterio SUM4.

### Precaución sobre `σ/√2`

En el centro, la interpretación como resolución equivalente de un extremo es razonable. En el borde, el extremo cercano y el lejano no tienen la misma estadística, así que el valor debe verse como una figura de mérito convencional del sistema L/R.

---

## Diapositiva 12 — F5: resolución SUM4 de EJ-230

La lectura es idéntica a la diapositiva 11.

### Valores

- `x = 0`: 139 ± 4 ps;
- `x = +400 mm`: 247 ± 7 ps;
- `x = +690 mm`: 419 ± 14 ps.

### Diferencia cualitativa

En el centro, EJ-230 y EJ-204 son prácticamente equivalentes dentro de las incertidumbres. Al avanzar hacia el borde, EJ-230 desarrolla:

- mayor anchura;
- cola más importante;
- menor eficiencia.

La eficiencia en `x = +690 mm` baja a 78.5 %, mucho más que el 95.5 % de EJ-204.

### Interpretación

La diferencia del borde es coherente con una menor robustez de la información del extremo lejano en la configuración simulada. La figura por sí sola no separa cuánto proviene de producción de luz, espectro de emisión, transporte, PDE o detalles del algoritmo SUM4.

---

## Diapositiva 13 — Comparación de resolución

### Lectura de la tabla

| Posición | EJ-204 | EJ-230 | Lectura |
|---|---:|---:|---|
| `x = 0` | 141 ± 4 ps | 139 ± 4 ps | equivalentes |
| `x = +400 mm` | 206 ± 6 ps | 247 ± 7 ps | EJ-230 más ancho |
| `x = +690 mm` | 365 ± 10 ps | 419 ± 14 ps | EJ-230 más ancho y menos eficiente |

### Porcentajes

En el borde:

- EJ-230 es aproximadamente `(419−365)/365 ≈ 14.8 %` mayor que EJ-204;
- dicho al revés, EJ-204 es aproximadamente `(419−365)/419 ≈ 12.9 %` menor que EJ-230.

La frase “15 % mejor” depende del denominador y conviene explicitarlo.

### Degradación centro-borde

- EJ-204: `(365−141)/141 ≈ 159 %`;
- EJ-230: `(419−139)/139 ≈ 201 %`.

### Significancia aproximada

Si las incertidumbres son independientes y puramente estadísticas:

- diferencia en `+400 mm`: `41 ± 9 ps`, cerca de 4.5 desviaciones estándar;
- diferencia en `+690 mm`: `54 ± 17 ps`, cerca de 3.1 desviaciones estándar.

Esto no incluye incertidumbres sistemáticas de modelado, topología SUM4 o electrónica.

### Mensaje correcto

Son resultados de Etapa 1. No son todavía la resolución final del detector real.

---

## Diapositiva 14 — F6: correlación entre SiPM TOP

### Pregunta que responde

¿Dos sensores superiores fluctúan juntos o aportan información complementaria?

### Ejes

Cada punto representa un evento:

- eje horizontal: `Npe` en el primer SiPM del par;
- eje vertical: `Npe` en el segundo;
- intensidad de color: número de eventos en cada bin, en escala logarítmica.

### Pearson `r`

- `r ≈ +1`: relación lineal positiva fuerte;
- `r ≈ 0`: poca relación lineal;
- `r < 0`: cuando uno aumenta, el otro tiende a disminuir.

Pearson mide solamente correlación lineal y puede ocultar estructuras de varias ramas.

### Par A: IDs 49/52

- posiciones aproximadas: −32 mm y +32 mm;
- separación: 64 mm;
- `r = −0.37`.

La anticorrelación indica comportamiento complementario a nivel agregado, pero no prueba por sí sola un mecanismo único. Las ramas visibles sugieren mezcla de geometrías o posiciones.

### Par B: IDs 50/51

- posiciones: −12 mm y +12 mm;
- cruzan el gap central de 24 mm;
- `r = +0.88`.

Ambos sensores responden de forma muy similar en el conjunto analizado. Son candidatos a información redundante.

### Par C: IDs 47/49

- mismo lado;
- separación: 40 mm;
- `r = +0.90`.

Sirve como control de correlación local esperada.

### Precaución importante

Una correlación alta no demuestra que se pueda eliminar un canal sin perder rendimiento. Para afirmar redundancia de diseño se necesitaría estudiar, por ejemplo:

- eficiencia de trigger al retirar un canal;
- resolución espacial;
- correlación condicionada a una misma posición `x_gun`;
- información mutua o desempeño de reconstrucción.

Si la correlación se calcula mezclando muchas posiciones, parte de `r` puede estar impulsada por la variación común con la posición y no por redundancia evento a evento.

---

## Diapositiva 15 — F7: estadística de orden en TOP

### Pregunta que responde

¿Cómo evoluciona el tiempo del primer, segundo, tercer, … fotón detectado en cada canal TOP y cuánto fluctúa entre eventos?

### Qué representa `n`

Dentro de cada evento y canal se ordenan los tiempos:

\[
t_{(1)}\le t_{(2)}\le t_{(3)}\le\cdots.
\]

`t_n` es el tiempo del fotón de orden `n`.

### Eje horizontal

`n = 1, 2, …, 20`. A medida que aumenta `n`, se seleccionan fotones progresivamente más tardíos, por lo que la curva media debe crecer.

### Eje vertical

\[
\langle t_n\rangle,
\]

promedio del tiempo del fotón `n` sobre todos los eventos válidos.

### Barras de error

No son el error de la media. Son la RMS evento a evento de `t_n`. Miden el ancho físico/estadístico de la distribución.

### Etiqueta `Npe`

Cada panel muestra el `⟨Npe⟩` medio de ese canal. Los canales más alejados de la zona de impacto suelen tener menos luz, tiempos medios más tardíos y mayor RMS.

### Piso estadístico

La aproximación usada es

\[
\sigma_{\mathrm{stat}}(t_n)
\approx
\frac{\sqrt{n}\,\tau_d}{\langle N_{pe}\rangle}.
\]

Es un límite idealizado para estadísticas de orden de una distribución exponencial, aproximadamente válido cuando `n` es pequeño frente a `Npe` y no existen dispersión geométrica ni otros términos.

### Cómo interpretar la banda

Si la RMS observada es mayor que `σ_stat`, existe variabilidad adicional que el modelo estadístico ideal no explica. Una candidata natural es la dispersión de caminos ópticos.

No se deben restar anchos linealmente. Si las contribuciones fueran independientes, la estimación correcta sería

\[
\sigma_{\mathrm{extra}}
\approx
\sqrt{
\sigma_{\mathrm{obs}}^2-
\sigma_{\mathrm{stat}}^2
}.
\]

Por tanto, la frase correcta es “el exceso sobre el piso evidencia una contribución adicional”, no “RMS menos piso es directamente la dispersión”.

### Geometría de los paneles

En `x = −690 mm` se muestran IDs 16–24, con posiciones de sensor desde aproximadamente −692 a −532 mm. El conjunto queda truncado por el borde izquierdo: no existen canales superiores más allá de ID 16.

### Evolución espacial

Al desplazarse de ID 16 hacia ID 24:

- disminuye en general `Npe`;
- aumenta el tiempo medio;
- aumenta la dispersión.

Esto es consistente con que los sensores más distantes reciben una mezcla de caminos más largos.

### Relación con F2

La relación es cualitativa:

- F7, en TOP, demuestra que existe dispersión adicional;
- F2, en END, no resuelve un piso constante `b` dentro del rango estudiado.

No son la misma población de sensores ni una medición directa del mismo estimador.

---

## Diapositiva 16 — Resumen y preguntas abiertas

### Resultados que sí están respaldados

1. En el centro, EJ-204 y EJ-230 tienen resolución equivalente dentro de las incertidumbres.
2. En posiciones alejadas del centro, EJ-230 degrada más y pierde más eficiencia.
3. La población que comienza en aproximadamente 4.7 ns en el borde está asociada al extremo lejano y es consistente con ToF geométrico.
4. El onset selecciona caminos más rápidos que el promedio y produce `v_eff,onset ≈ 292 mm/ns`.
5. La resolución se correlaciona mejor con la luz del extremo débil que con la luz promedio de ambos extremos.
6. F7 muestra una RMS superior al piso ideal de estadística de orden, compatible con dispersión adicional de caminos.
7. La fila TOP mantiene una señal muy por encima de 20 PE en la simulación EJ-204 EndTop.

### Formulaciones que conviene suavizar

- En vez de “Poisson puro”, decir “tendencia dominada por estadística de fotones, sin piso constante resuelto en este rango”.
- En vez de “RMS − floor es dispersión”, decir “el exceso cuadrático sobre el piso es compatible con dispersión adicional”.
- En vez de “los sensores correlacionados son prescindibles”, decir “son candidatos a redundancia y requieren una prueba de desempeño”.
- En F3, no llamar “simétrico” a un evento central con 480 vs 123 hits.

### Pendientes reales

- definición exacta de SUM4;
- significado definitivo de T4/T20;
- tratamiento del jitter de 20 ps;
- inclusión de SPTR, time-walk y electrónica;
- análisis EndTop de EJ-230;
- separación de reflexión y emisión lenta en la cola central de EJ-230;
- correlaciones TOP condicionadas a posición.

---

# 4. Hilo argumental recomendado para exponer

Una secuencia clara para presentar es:

1. **Geometría:** definir END_LEFT, END_RIGHT y TOP.
2. **F5:** mostrar el hecho principal: la resolución empeora hacia el borde.
3. **F3:** visualizar que en el borde casi toda la luz llega al extremo cercano.
4. **F2:** demostrar que la variable útil es la luz del extremo débil, no el promedio.
5. **F1 + QA-1c:** explicar que la estructura tardía del borde es ToF del extremo lejano.
6. **F7:** mostrar que existe dispersión óptica adicional, aunque F2 no resuelva un piso constante.
7. **F4/F6:** discutir la cobertura y posible redundancia del sistema TOP.
8. **Comparación final:** EJ-204 y EJ-230 son similares en el centro, pero EJ-204 es más robusto al alejarse del centro en esta simulación.

Una frase de cierre técnicamente prudente sería:

> “En esta Etapa 1, el deterioro del timing hacia los bordes está fuertemente asociado a la escasa estadística del extremo lejano. La dispersión de caminos existe, pero el ajuste simple de F2 no resuelve todavía un piso constante distinto de cero. La comparación debe completarse con la cadena electrónica y con una definición validada del estimador SUM4.”
