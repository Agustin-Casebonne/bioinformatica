# Práctica 2 (Sesión 3): alineamientos con PairwiseAligner

Ejecutar desde esta carpeta, con el entorno virtual activado:

```bash
pip install biopython pandas matplotlib seaborn
python ejercicio1_grid_search.py
python ejercicio2_heatmap.py
python ejercicio3_tiempos.py
```

## Ejercicio 1: Grid Search

**Espacio de búsqueda.** `open_gap` de -1 a -20 con paso 1 y `extend_gap` de -0.1 a -5.0 con paso 0.1, es decir, 1000 combinaciones. Se usa alineamiento global, con match = 1 y mismatch = 0.

**Resultados.**

| Relación entre parámetros | Combinaciones | Empates |
|---|---|---|
| `open_gap < extend_gap` (abrir es más caro) | 895 | siempre 1 |
| `open_gap == extend_gap` | 5 | 9 |
| `open_gap > extend_gap` (abrir es más barato) | 100 | de 5 a 111 |

**Explicación matemática.** `seq_sana` tiene 14 bases y `seq_paciente` 8, así que hay que colocar al menos 6 posiciones de hueco. Si las 8 bases coinciden y esas 6 posiciones se reparten en `k` huecos, la puntuación es:

```
score = 8 + k * open_gap + (6 - k) * extend_gap
```

Cada hueco adicional cambia la puntuación en `open_gap - extend_gap`:

- Si esa diferencia es negativa, conviene tener un solo hueco (`k = 1`). Solo hay una forma de hacerlo con 8 coincidencias, y por eso la solución es única.
- Si es cero, `k` no afecta a la puntuación y aparecen los empates.
- Si es positiva, compensa fragmentar el hueco, y hay muchas formas de hacerlo.

**Punto de colapso:** la frontera exacta es `open_gap = extend_gap`. La solución es única justo cuando `open_gap < extend_gap`, es decir, cuando abrir un hueco penaliza más que extenderlo.

**Solución tras el colapso** (`open_gap = -2`, `extend_gap = -1`, score = 8 - 2 - 5 = 1):

```
target  GCATGCATGCATGC
        ||||||------||
query   GCATGC------GC
```

**Aportación adicional.** `ejercicio1_heatmap.png` muestra el número de empates en todo el espacio. Se ve la escalera que separa las dos regiones. El máximo (111 empates) aparece en `extend_gap = open_gap - 1`: ahí, cada hueco extra gana 1 punto, lo mismo que vale una coincidencia.

**Lectura biológica.** Tener una solución única no significa que sea la correcta. Solo indica que los parámetros deshacen el empate. La pérdida de un único bloque de 6 bases es una hipótesis razonable, porque una deleción grande suele ser un solo evento, pero el alineamiento no lo demuestra. Además, esta secuencia es repetitiva (`GCAT`), lo que facilita los empates. Por eso los parámetros no deben elegirse solo para que salga un resultado único: hacerlo sería el peligro que señala la nota del enunciado.

## Ejercicio 2: comparación de 5 especies

**Datos.** `citocromo_c.fasta` contiene el citocromo c descargado de UniProt (`https://rest.uniprot.org/uniprotkb/<ID>.fasta`):

| Especie | Identificador |
|---|---|
| Humano | P99999 |
| Chimpancé | P99998 |
| Ratón | P62897 |
| Pez cebra | Q6IQM2 |
| Mosca (*Drosophila*) | P04657 |

**Parámetros.** Alineamiento global, con `open_gap = -10` y `extend_gap = -0.5`.

**Resultados.** Score con BLOSUM62 (`ejercicio2_heatmap_BLOSUM62.png`):

| | Humano | Chimpancé | Ratón | Pez cebra | Mosca |
|---|---|---|---|---|---|
| **Humano** | 565 | 565 | 522 | 475 | 368 |
| **Chimpancé** | 565 | 565 | 522 | 475 | 368 |
| **Ratón** | 522 | 522 | 565 | 504 | 389 |
| **Pez cebra** | 475 | 475 | 504 | 558 | 377.5 |
| **Mosca** | 368 | 368 | 389 | 377.5 | 565 |

**Comentarios.**

- La matriz es simétrica y el valor más alto de cada fila está en la diagonal, como cabe esperar.
- Humano y chimpancé tienen scores idénticos porque su citocromo c es la misma secuencia.
- Desde el humano, el score baja de ratón a pez y a mosca, en línea con la distancia evolutiva.
- El score no es un reloj evolutivo exacto: humano y ratón están igual de lejos del pez, pero puntúan distinto (475 y 504).

**Otras matrices.** Parecido con el humano, calculado como score / score(Humano-Humano):

| Matriz | Ratón | Pez cebra | Mosca |
|---|---|---|---|
| BLOSUM45 | 0.92 | 0.84 | 0.67 |
| BLOSUM62 | 0.92 | 0.84 | 0.65 |
| BLOSUM80 | 0.92 | 0.83 | 0.64 |
| PAM250 | 0.94 | 0.87 | 0.70 |
| PAM30 | 0.89 | 0.79 | 0.51 |

- Los scores absolutos no se pueden comparar entre matrices porque cada una tiene su escala. Por eso se normalizan.
- El orden de las especies no cambia con ninguna matriz.
- Las matrices pensadas para secuencias lejanas (PAM250, BLOSUM45) son más tolerantes con las sustituciones y dan valores relativos más altos.
- Las pensadas para secuencias cercanas (PAM30, BLOSUM80) son más estrictas. PAM30 es la que más separa a las especies.

## Ejercicio 3: escalabilidad

**Método.** Secuencias aleatorias generadas con `random` y tiempo medido con `time.perf_counter`. Cada tamaño se mide 3 veces y se guarda el mejor tiempo, para reducir el ruido.

**Resultados.** Tiempos medidos en esta ejecución (varían según el ordenador):

| Longitud (bases) | 100 | 500 | 1000 | 2000 | 5000 |
|---|---|---|---|---|---|
| Tiempo (ms) | 0.03 | 1.26 | 5.19 | 22.51 | 142.21 |

**Interpretación.** `ejercicio3_tiempos.png` muestra una curva, no una recta. Al duplicar la longitud de 1000 a 2000, el tiempo se multiplica por 4.3, y de 2000 a 5000 (x 2.5) se multiplica por 6.3. Ambos valores están muy cerca de lo esperado (4 y 6.25), lo que confirma un coste **O(n²)**: la matriz tiene n x m celdas y hay que rellenarlas todas. El salto de 100 a 500 (x 47 en lugar de x 25) se aleja más, pero con 0.03 ms la medida es muy sensible al ruido.

**¿Escala hasta 1 millón de bases?** No para un pipeline con muchos alineamientos:

- Con 1 millón de bases, la matriz tendría 10¹² celdas.
- Extrapolando desde 5000 bases, cada par tardaría unas **1.6 horas**.
- Es solo una estimación, y no tiene en cuenta la memoria, que también puede ser un límite con una matriz de ese tamaño.

Por eso, a gran escala se usan métodos heurísticos como BLAST, que renuncian a garantizar el óptimo a cambio de velocidad.
