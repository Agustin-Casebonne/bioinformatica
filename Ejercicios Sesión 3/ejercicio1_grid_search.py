"""Ejercicio 1: Grid Search de las penalizaciones de hueco de PairwiseAligner."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import Align

carpeta = Path(__file__).resolve().parent

seq_sana = "GCATGCATGCATGC"
seq_paciente = "GCATGCGC"

# Espacio de búsqueda elegido:
#   open_gap:   de -1 a -20, de 1 en 1      (20 valores)
#   extend_gap: de -0.1 a -5.0, de 0.1 en 0.1 (50 valores)
valores_open = [-float(i) for i in range(1, 21)]
valores_extend = [-i / 10 for i in range(1, 51)]


def crear_alineador(open_gap, extend_gap):
    alineador = Align.PairwiseAligner()
    alineador.mode = "global"
    alineador.match_score = 1.0
    alineador.mismatch_score = 0.0
    alineador.open_gap_score = open_gap
    alineador.extend_gap_score = extend_gap
    return alineador


# 1. Grid Search con dos bucles anidados
filas = []
for open_gap in valores_open:
    for extend_gap in valores_extend:
        alineador = crear_alineador(open_gap, extend_gap)
        alineamientos = alineador.align(seq_sana, seq_paciente)
        filas.append({
            "open_gap": open_gap,
            "extend_gap": extend_gap,
            "empates": len(alineamientos),
            "score": alineamientos.score,
        })


df = pd.DataFrame(filas)
df.to_csv(carpeta / "ejercicio1_resultados.csv", index=False)
print("Combinaciones probadas:", len(df))
print(df.head(), "\n")

# 2. Punto de colapso: comparamos el número de empates según la relación
#    entre abrir un hueco y extenderlo
df["relacion"] = np.select(
    [df["open_gap"] > df["extend_gap"], df["open_gap"] == df["extend_gap"]],
    ["abrir MENOS caro que extender", "abrir IGUAL que extender"],
    default="abrir MÁS caro que extender",
)
resumen = df.groupby("relacion")["empates"].agg(["count", "min", "max"])
print("Empates según la relación entre open_gap y extend_gap:")
print(resumen, "\n")

# Para cada extend_gap, el primer open_gap (bajando desde -1) con solución única
print("Frontera: primer open_gap con solución única para algunos extend_gap")
for extend_gap in [-0.5, -1.0, -2.0, -3.0, -4.0, -5.0]:
    tramo = df[df["extend_gap"] == extend_gap]
    print(f"  extend_gap={extend_gap:5}  ->  open_gap="
          f"{tramo[tramo['empates'] == 1]['open_gap'].max()}")

# 3. Solución obtenida tras el colapso
alineador = crear_alineador(-2.0, -1.0)
alineamientos = alineador.align(seq_sana, seq_paciente)
print("\nSolución tras el colapso (open_gap=-2, extend_gap=-1):")
print("Empates:", len(alineamientos), " Score:", alineamientos.score)
print(alineamientos[0])

# Mismo problema justo en la frontera, para ver los empates
alineador = crear_alineador(-1.0, -1.0)
alineamientos = alineador.align(seq_sana, seq_paciente)
print("En la frontera (open_gap=-1, extend_gap=-1):")
print("Empates:", len(alineamientos), " Score:", alineamientos.score)
print("Dos ejemplos de empate:")
print(alineamientos[0])
print(alineamientos[1])

# 4. Aportación adicional: mapa de calor de todo el espacio de parámetros
tabla = df.pivot(index="open_gap", columns="extend_gap", values="empates")
tabla = tabla.sort_index(ascending=False).sort_index(axis=1, ascending=False)

plt.figure(figsize=(12, 6))
sns.heatmap(np.log10(tabla), cmap="viridis", xticklabels=5,
            cbar_kws={"label": "log10(número de empates)"})
plt.title("Empates óptimos según las penalizaciones de hueco")
plt.xlabel("extend_gap")
plt.ylabel("open_gap")
plt.tight_layout()
plt.savefig(carpeta / "ejercicio1_heatmap.png", dpi=150)
print("Gráfico guardado en ejercicio1_heatmap.png")
