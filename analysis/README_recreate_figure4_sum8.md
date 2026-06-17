# Macro ROOT: recreación de la Figura 4 y velocidad efectiva

Archivo principal: `recreate_figure4_sum8.C`

## Qué analiza

La macro procesa el `TTree sipm_hits` de cuatro campañas:

| Campaña | Directorio predeterminado |
|---|---|
| EJ-204, sólo END | `/home/reriosto/SHiP/t0minidaq/endonly_mylar_20260614` |
| EJ-230, sólo END | `/home/reriosto/SHiP/t0minidaq/endonly_mylar_230` |
| EJ-204, END+TOP | `/home/reriosto/SHiP/t0minidaq/sslg4/exec07_endtop_2000` |
| EJ-230, END+TOP | `/home/reriosto/SHiP/t0minidaq/results_ej230/data` |

Busca automáticamente los archivos `photon_hits*.root` dentro de cada directorio y toma la posición desde la rama `gun_x_mm`; por ello admite tanto nombres del tipo `photon_hits_x400mm.root` como nombres por número de run.

## Definiciones de lectura

- **END izquierdo:** suma analógica de los ocho SiPM con `global_id = 0...7`.
- **END derecho:** suma analógica de los ocho SiPM con `global_id = 8...15`.
- **TOP8:** suma analógica de los ocho SiPM TOP geométricamente más próximos a `gun_x_mm`.
- **TOP-all:** suma de los 70 SiPM TOP (`global_id = 16...85`), guardada como diagnóstico adicional.

TOP8 es la curva principal para TOP porque la geometría superior es distribuida a lo largo de la barra y no constituye un único arreglo físico de ocho sensores como el END. El número de sensores vecinos es configurable.

Cada fotoelectrón aporta una respuesta de un PE modelada mediante una diferencia de exponenciales normalizada, con `tau_r = 0.5 ns` y `tau_f = 5 ns`. El tiempo corresponde al primer cruce de un umbral leading-edge de `4 PE` por defecto.

## Ejecución

Copiar la macro, por ejemplo, en:

```bash
cp recreate_figure4_sum8.C /home/reriosto/SHiP/orchestrator/analysis/
cd /home/reriosto/SHiP/orchestrator
```

Ejecutar con ACLiC:

```bash
root -l -b -q 'analysis/recreate_figure4_sum8.C+()'
```

La salida predeterminada es:

```text
/home/reriosto/SHiP/orchestrator/outputs/figure4_sum8_validation.root
/home/reriosto/SHiP/orchestrator/outputs/figure4_sum8_validation.csv
/home/reriosto/SHiP/orchestrator/outputs/figure4_sum8_validation_plots/
```

Para cambiar rutas o parámetros:

```bash
root -l -b -q 'analysis/recreate_figure4_sum8.C+("/ruta/salida.root", \
  "/ruta/ej204_endonly", "/ruta/ej230_endonly", \
  "/ruta/ej204_endtop", "/ruta/ej230_endtop", \
  4.0, 8, 2000)'
```

Los últimos parámetros son:

1. umbral leading-edge en PE;
2. número de SiPM TOP vecinos que se suman;
3. número esperado de eventos por posición.

## Contenido del archivo ROOT

Para cada dataset hay un directorio con:

- `summary`: `TTree` con medias, resoluciones, eficiencias y canales TOP seleccionados;
- `time_distributions/`: histogramas `TH1D` por posición;
- gráficos `TGraphErrors` de resolución, tiempo medio y eficiencia;
- ajustes lineales `TF1` usados para la velocidad efectiva;
- `TCanvas` completos:
  - `c_fig4_*`: equivalente de la Figura 4 para END SUM8;
  - `c_top_*`: resolución TOP8 y TOP-all;
  - `c_speed_*`: tiempos medios y cálculo de velocidad efectiva;
  - `c_efficiency_*`: eficiencia de cruce del umbral.

En la raíz del archivo se guarda además `velocity_summary`, con los resultados de todos los datasets.

## Cálculo de la velocidad efectiva

La macro ajusta tres relaciones lineales usando `x` en centímetros:

```text
<t_L>       = a_L + b_L x       => v_L     = 1/|b_L|
<t_R>       = a_R + b_R x       => v_R     = 1/|b_R|
<t_L-t_R>   = a_D + b_D x       => v_delta = 2/|b_D|
```

El resultado `v_eff_consensus_cm_ns` es una combinación ponderada **diagnóstica** de los tres estimadores cuando sus incertidumbres son válidas. No debe interpretarse como una combinación estadística independiente, porque los tres estimadores usan los mismos eventos. Para la validación principal conviene citar `v_delta`, y usar `v_L` y `v_R` como verificaciones cruzadas. En los canvases se compara con `15.5 cm/ns`.

## Observaciones para la validación

- La velocidad debe extraerse principalmente de los tiempos END; TOP8 es una lectura local y su tiempo no tiene por qué presentar la misma pendiente a lo largo de toda la barra.
- La macro usa directamente `time_ns`, que en estas simulaciones contiene tiempo de emisión, propagación óptica y el jitter por hit configurado en Geant4.
- La resolución END individual se obtiene directamente de `sigma(t_L)` y `sigma(t_R)`, porque el tiempo primario de Geant4 sirve como referencia común. También se guardan la resolución de `(t_L+t_R)/2`, el promedio ponderado y `t_L-t_R`.
- El promedio ponderado usa pesos `1/sigma^2` por posición y no incluye explícitamente una covarianza entre ambos extremos.
