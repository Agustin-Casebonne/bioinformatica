"""Ejercicio 0: Needleman-Wunsch y Smith-Waterman en Python puro.

a) Los dos algoritmos usan solo Python estándar (listas), siguiendo
   los 4 pasos de la presentación:
   1. Puntuaciones  2. Iniciación de la matriz  3. Rellenado  4. Rastreo inverso
b) Se comparan los resultados con PairwiseAligner de Biopython.
"""
import random

from Bio.Align import PairwiseAligner


# =========================================================
# a) ALGORITMOS EN PYTHON PURO
# =========================================================

def puntuacion(x, y, match, mismatch):
    """Paso 1: S(x, y) vale 'match' si las letras coinciden y 'mismatch' si no."""
    if x == y:
        return match
    return mismatch


def alinear(seq1, seq2, match, mismatch, hueco, modo):
    """Alinea seq1 (filas) con seq2 (columnas). modo = 'global' o 'local'."""
    n = len(seq1)
    m = len(seq2)

    # Paso 2. Iniciación de la matriz (n+1) x (m+1)
    C = [[0] * (m + 1) for _ in range(n + 1)]
    if modo == "global":
        # Needleman-Wunsch: primera fila y columna acumulan la penalización de hueco
        for i in range(1, n + 1):
            C[i][0] = i * hueco
        for j in range(1, m + 1):
            C[0][j] = j * hueco
    # Smith-Waterman: primera fila y columna se quedan a cero

    # Paso 3. Rellenado de la matriz
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            v1 = C[i - 1][j - 1] + puntuacion(seq1[i - 1], seq2[j - 1], match, mismatch)  # diagonal
            v2 = C[i - 1][j] + hueco  # desde arriba
            v3 = C[i][j - 1] + hueco  # desde la izquierda
            v4 = max(v1, v2, v3)
            if modo == "local" and v4 < 0:
                v4 = 0  # Smith-Waterman trunca los negativos a cero
            C[i][j] = v4

    # Paso 4. Rastreo inverso: elegir el punto de inicio
    if modo == "global":
        # Needleman-Wunsch: esquina inferior derecha
        i, j = n, m
    else:
        # Smith-Waterman: celda con el valor más alto
        i, j = 0, 0
        for fila in range(n + 1):
            for col in range(m + 1):
                if C[fila][col] > C[i][j]:
                    i, j = fila, col
    score = C[i][j]

    # Reconstruimos el alineamiento de atrás hacia delante
    alin1 = ""
    alin2 = ""
    while i > 0 or j > 0:
        if modo == "local" and C[i][j] == 0:
            break  # Smith-Waterman termina al llegar a un cero

        if i > 0 and j > 0 and C[i][j] == C[i - 1][j - 1] + puntuacion(seq1[i - 1], seq2[j - 1], match, mismatch):
            # Diagonal: una letra de cada secuencia
            alin1 = seq1[i - 1] + alin1
            alin2 = seq2[j - 1] + alin2
            i, j = i - 1, j - 1
        elif i > 0 and C[i][j] == C[i - 1][j] + hueco:
            # Arriba: letra de la secuencia 1 y hueco en la secuencia 2
            alin1 = seq1[i - 1] + alin1
            alin2 = "-" + alin2
            i = i - 1
        else:
            # Izquierda: letra de la secuencia 2 y hueco en la secuencia 1
            alin1 = "-" + alin1
            alin2 = seq2[j - 1] + alin2
            j = j - 1

    return C, score, alin1, alin2


def needleman_wunsch(seq1, seq2, match, mismatch, hueco):
    return alinear(seq1, seq2, match, mismatch, hueco, "global")


def smith_waterman(seq1, seq2, match, mismatch, hueco):
    return alinear(seq1, seq2, match, mismatch, hueco, "local")


# =========================================================
# FUNCIONES PARA MOSTRAR Y COMPARAR RESULTADOS
# =========================================================

def imprimir_matriz(C, seq1, seq2):
    print("       " + "".join(f"{letra:>5}" for letra in seq2))
    for i, fila in enumerate(C):
        letra = " " if i == 0 else seq1[i - 1]
        print(f"{letra:>2}" + "".join(f"{valor:>5}" for valor in fila))


def linea_simbolos(alin1, alin2):
    """'|' coincidencia, '*' desemparejamiento, '-' hueco (como en la presentación)."""
    simbolos = ""
    for a, b in zip(alin1, alin2):
        if a == "-" or b == "-":
            simbolos += "-"
        elif a == b:
            simbolos += "|"
        else:
            simbolos += "*"
    return simbolos


def comparar_con_biopython(seq1, seq2, match, mismatch, hueco, modo):
    print("=" * 50)
    print(f" MODO {modo.upper()}: {seq1} vs {seq2}")
    print(f" match={match}, mismatch={mismatch}, hueco={hueco}")
    print("=" * 50)

    # 1. Nuestra implementación
    C, score, alin1, alin2 = alinear(seq1, seq2, match, mismatch, hueco, modo)
    print("Matriz de programación dinámica:")
    imprimir_matriz(C, seq1, seq2)
    print(f"\n>>> NUESTRA IMPLEMENTACIÓN")
    print(f"Puntuación: {score}")
    print(alin1)
    print(linea_simbolos(alin1, alin2))
    print(alin2)

    # 2. Biopython con la misma configuración
    alineador = PairwiseAligner()
    alineador.mode = modo
    alineador.match_score = match
    alineador.mismatch_score = mismatch
    alineador.open_gap_score = hueco
    alineador.extend_gap_score = hueco  # penalización lineal, como en la presentación

    alineamientos = alineador.align(seq1, seq2)
    print(f"\n>>> BIOPYTHON")
    print(f"Puntuación: {alineamientos.score}")
    print(f"Número de alineamientos óptimos: {len(alineamientos)}")
    print(alineamientos[0])

    # 3. Comparación
    optimos_bio = [(a[0], a[1]) for a in alineamientos]
    print(">>> COMPARACIÓN")
    print(f"¿Misma puntuación?: {score == alineamientos.score}")
    print(f"¿Nuestro alineamiento es uno de los óptimos de Biopython?: "
          f"{(alin1, alin2) in optimos_bio}\n")


# =========================================================
# b) PRUEBAS
# =========================================================

# Prueba 1: ejemplo de la presentación (debe dar las mismas matrices que las diapositivas)
comparar_con_biopython("AGCGTTA", "ACGTGA", 5, -4, -6, "global")
comparar_con_biopython("CGT", "ACGTGA", 5, -4, -6, "local")

# Prueba 2: ejemplo local del notebook de clase (la isla ACG dentro de GACGG)
comparar_con_biopython("GACGG", "ACG", 2, -1, -1, "global")
comparar_con_biopython("GACGG", "ACG", 2, -1, -1, "local")

# Prueba 3: secuencia con un empate (pérdida de una T)
comparar_con_biopython("GATTACA", "GATACA", 5, -4, -6, "global")

# Prueba 4: aminoácidos generados aleatoriamente (longitud 10 o menos)
random.seed(1)
aminoacidos = "ACDEFGHIKLMNPQRSTVWY"
prot1 = "".join(random.choices(aminoacidos, k=8))
prot2 = "".join(random.choices(aminoacidos, k=6))
comparar_con_biopython(prot1, prot2, 2, -1, -2, "global")
comparar_con_biopython(prot1, prot2, 2, -1, -2, "local")
