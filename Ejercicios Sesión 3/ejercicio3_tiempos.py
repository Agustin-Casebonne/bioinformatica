"""Ejercicio 3: cuánto tarda PairwiseAligner según la longitud de las secuencias."""
import random
import time
from pathlib import Path

import matplotlib.pyplot as plt
from Bio import Align

carpeta = Path(__file__).resolve().parent

random.seed(42)  # para que las secuencias sean siempre las mismas
tamanos = [100, 500, 1000, 2000, 5000]
repeticiones = 3  # se repite cada medida y se guarda la mejor, para reducir el ruido


def secuencia_aleatoria(longitud):
    return "".join(random.choices("ACGT", k=longitud))


alineador = Align.PairwiseAligner()
alineador.mode = "global"

tiempos = []
for n in tamanos:
    seq1 = secuencia_aleatoria(n)
    seq2 = secuencia_aleatoria(n)

    medidas = []
    for _ in range(repeticiones):
        inicio = time.perf_counter()
        alineamientos = alineador.align(seq1, seq2)
        primero = alineamientos[0]  # fuerza a reconstruir un alineamiento óptimo
        fin = time.perf_counter()
        medidas.append((fin - inicio) * 1000)  # de segundos a milisegundos

    tiempos.append(min(medidas))
    print(f"n = {n:5}  ->  {tiempos[-1]:10.2f} ms")

# Si el coste es O(n²), al multiplicar n por k el tiempo se multiplica por k²
print("\nComprobación del crecimiento:")
for i in range(1, len(tamanos)):
    k = tamanos[i] / tamanos[i - 1]
    print(f"n x {k:.1f}: tiempo x {tiempos[i] / tiempos[i - 1]:6.1f} "
          f"(esperado si O(n²): {k ** 2:5.1f})")

# Curva teórica n² ajustada al punto más grande
constante = tiempos[-1] / tamanos[-1] ** 2
teorica = [constante * n ** 2 for n in tamanos]

# Extrapolación a 1 millón de bases
horas = constante * 1_000_000 ** 2 / 1000 / 3600
print(f"\nEstimación para n = 1.000.000: {horas:.1f} horas")
print(f"Celdas de la matriz con n = 1.000.000: {1_000_000 ** 2:.0e}")

plt.figure(figsize=(8, 5))
plt.plot(tamanos, tiempos, "o-", label="Tiempo medido")
plt.plot(tamanos, teorica, "--", label="Curva teórica n²")
plt.xlabel("Longitud de la secuencia (bases)")
plt.ylabel("Tiempo de ejecución (ms)")
plt.title("Escalabilidad de PairwiseAligner")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(carpeta / "ejercicio3_tiempos.png", dpi=150)
print("Gráfico guardado en ejercicio3_tiempos.png")
