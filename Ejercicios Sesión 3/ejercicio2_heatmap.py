"""Ejercicio 2: comparación "todos contra todos" del citocromo c en 5 especies."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import Align, SeqIO
from Bio.Align import substitution_matrices

carpeta = Path(__file__).resolve().parent

# Mismo orden en el que están las secuencias en citocromo_c.fasta
especies = ["Humano", "Chimpancé", "Ratón", "Pez cebra", "Mosca"]

# 1 Cargar el FASTA con Bio.SeqIO
registros = list(SeqIO.parse(carpeta / "citocromo_c.fasta", "fasta"))
for especie, registro in zip(especies, registros):
    print(f"{especie:10} {registro.id:25} {len(registro.seq)} aminoácidos")


def matriz_scores(nombre_matriz):
    """Devuelve un DataFrame con el score global de cada par de especies."""
    alineador = Align.PairwiseAligner()
    alineador.mode = "global"
    alineador.substitution_matrix = substitution_matrices.load(nombre_matriz)
    alineador.open_gap_score = -10
    alineador.extend_gap_score = -0.5

    # 2 Bucle de todos contra todos: 5 x 5 = 25 puntuaciones
    scores = np.zeros((5, 5))
    for i in range(5):
        for j in range(5):
            scores[i, j] = alineador.score(registros[i].seq, registros[j].seq)

    # 3 Volcado a un DataFrame
    return pd.DataFrame(scores, index=especies, columns=especies)


# 4 Heatmap con BLOSUM62 y experimentos con otras matrices
matrices = ["BLOSUM62", "BLOSUM45", "BLOSUM80", "PAM250", "PAM30"]

for nombre in matrices:
    tabla = matriz_scores(nombre)
    print(f"\nScores con {nombre}:")
    print(tabla.round(1))

    plt.figure(figsize=(7, 5.5))
    sns.heatmap(tabla, annot=True, fmt=".0f", cmap="rocket",
                cbar_kws={"label": "Score"})
    plt.title(f"Citocromo c: score de alineamiento global ({nombre})")
    plt.tight_layout()
    plt.savefig(carpeta / f"ejercicio2_heatmap_{nombre}.png", dpi=150)
    plt.close()

    # Normalizamos por el score de la proteína humana contra sí misma
    # para poder comparar matrices con escalas distintas
    humano = tabla.loc["Humano"]
    print("Parecido relativo a Humano-Humano:")
    print((humano / humano["Humano"]).round(2).to_string())

print("\nGráficos guardados como ejercicio2_heatmap_<matriz>.png")
