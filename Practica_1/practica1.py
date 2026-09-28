import random

# 1. TABLAS DE REFERENCIA BIOLOGICA

COMPLEMENTO_ADN = {"A": "T", "T": "A", "C": "G", "G": "C"}

# Complementariedad entre ADN y ARN
COMPLEMENTO_ARN = {"A": "U", "T": "A", "C": "G", "G": "C"}

BASES_VALIDAS = set("ATCG")

CODIGO_GENETICO = {
    "UUU": "Fenilalanina", "UUC": "Fenilalanina",
    "UUA": "Leucina", "UUG": "Leucina",
    "CUU": "Leucina", "CUC": "Leucina", "CUA": "Leucina", "CUG": "Leucina",
    "AUU": "Isoleucina", "AUC": "Isoleucina", "AUA": "Isoleucina",
    "AUG": "Metionina",
    "GUU": "Valina", "GUC": "Valina", "GUA": "Valina", "GUG": "Valina",
    "UCU": "Serina", "UCC": "Serina", "UCA": "Serina", "UCG": "Serina",
    "CCU": "Prolina", "CCC": "Prolina", "CCA": "Prolina", "CCG": "Prolina",
    "ACU": "Treonina", "ACC": "Treonina", "ACA": "Treonina", "ACG": "Treonina",
    "GCU": "Alanina", "GCC": "Alanina", "GCA": "Alanina", "GCG": "Alanina",
    "UAU": "Tirosina", "UAC": "Tirosina",
    "CAU": "Histidina", "CAC": "Histidina",
    "CAA": "Glutamina", "CAG": "Glutamina",
    "AAU": "Asparagina", "AAC": "Asparagina",
    "AAA": "Lisina", "AAG": "Lisina",
    "GAU": "Acido aspartico", "GAC": "Acido aspartico",
    "GAA": "Acido glutamico", "GAG": "Acido glutamico",
    "UGU": "Cisteina", "UGC": "Cisteina",
    "UGG": "Triptofano",
    "CGU": "Arginina", "CGC": "Arginina", "CGA": "Arginina", "CGG": "Arginina",
    "AGU": "Serina", "AGC": "Serina", "AGA": "Arginina", "AGG": "Arginina",
    "GGU": "Glicina", "GGC": "Glicina", "GGA": "Glicina", "GGG": "Glicina",
}

# Codones que indican el final de la traducción
STOP_CODONES = {"UAA": "Ocre", "UAG": "Ambar", "UGA": "Opalo"}

CODON_INICIO = "AUG"


def linea(titulo=None):
    print()
    if titulo:
        print(titulo)
        print()


def complementaria_adn(cadena):
    """Devuelve la hebra complementaria de una cadena de ADN."""
    return "".join(COMPLEMENTO_ADN[b] for b in cadena)


def complementaria_inversa_adn(cadena):
    """Devuelve la complementaria inversa de una cadena de ADN."""
    return complementaria_adn(cadena)[::-1]


def es_adn_valido(cadena):
    return len(cadena) > 0 and all(b in BASES_VALIDAS for b in cadena)


def mostrar_doble_helice(hebra_5_3, hebra_3_5, titulo="Doble helice"):
    print(f"\n{titulo}")
    print("  5'-" + hebra_5_3 + "-3'")
    print("     " + "|" * len(hebra_5_3))
    print("  3'-" + hebra_3_5 + "-5'")


# 2. OBTENER LA SECUENCIA DE ADN INICIAL

def generar_cadena_aleatoria(num_codones=6):
    """
    Genera una secuencia de ADN que comienza con AUG y termina
    con un codón de parada.
    """

    codones_intermedios = [c for c in CODIGO_GENETICO if c != CODON_INICIO]
    elegidos = [random.choice(codones_intermedios) for _ in range(num_codones)]
    codon_stop = random.choice(list(STOP_CODONES.keys()))

    arn_referencia = CODON_INICIO + "".join(elegidos) + codon_stop

    # La hebra codificante de ADN coincide con el ARN, cambiando U por T
    hebra_codificante = arn_referencia.replace("U", "T")

    # Se obtiene la hebra complementaria
    hebra_molde = complementaria_adn(hebra_codificante)

    return hebra_codificante, hebra_molde


def pedir_cadena_manual():
    while True:
        entrada = input(
            "Introduce la hebra codificante de ADN (solo A,T,C,G): "
        ).strip().upper()

        if es_adn_valido(entrada):
            return entrada, complementaria_adn(entrada)

        print("  Secuencia no valida. Usa solo las letras A, T, C, G.")


def leer_cadena_de_fichero(ruta):
    """Lee una secuencia de ADN desde un fichero de texto o FASTA."""

    with open(ruta) as f:
        lineas_fichero = [
            l.strip()
            for l in f
            if l.strip() and not l.startswith(">")
        ]

    secuencia = "".join(lineas_fichero).upper()

    if not es_adn_valido(secuencia):
        raise ValueError(
            "El fichero no contiene una secuencia de ADN valida "
            "(solo A,T,C,G)."
        )

    return secuencia, complementaria_adn(secuencia)


# 3. REPLICACIÓN DEL ADN

TAM_CEBADOR = 3      # bases del cebador de ARN (en la realidad, unas 10)
TAM_FRAGMENTO = 6    # bases por fragmento de Okazaki (en la realidad, 100-2000)


def cebador_arn(fragmento_adn, largo):
    """Cebador de ARN: primeras bases de la nueva hebra, con U en lugar de T."""
    return fragmento_adn[:largo].replace("T", "U")


def replicacion(hebra_5_3, hebra_3_5):
    """
    Simula la replicacion del ADN.
    Se muestran la cadena lider, que se sintetiza de forma continua,
    y la cadena rezagada, que se sintetiza mediante fragmentos de Okazaki.
    Cada sintesis empieza con un cebador de ARN.
    """

    linea("PASO 1: REPLICACIÓN DEL ADN  (ADN -> ADN)")

    mostrar_doble_helice(hebra_5_3, hebra_3_5, "ADN original:")

    print("\n> HELICASA: separa las dos hebras del ADN (horquilla de replicacion).")

    # Hebra lider

    print("\n> HEBRA LIDER")
    print("  Se sintetiza de forma continua, con un unico cebador.")

    hebra_lider = complementaria_adn(hebra_3_5)
    k = min(TAM_CEBADOR, len(hebra_lider))
    cebador = cebador_arn(hebra_lider, k)

    print(f"  Molde (3'->5'): {hebra_3_5}")
    print(f"  PRIMASA: sintetiza el cebador de ARN  5'-{cebador}-3'")
    print("  ADN POLIMERASA III: alarga la hebra desde el cebador en sentido 5'->3'.")
    print(f"  Hebra en sintesis (5'->3'): {cebador}{hebra_lider[k:]}  (ARN + ADN)")
    print("  ADN POLIMERASA I: sustituye el cebador de ARN por ADN.")
    print(f"  Nueva hebra (5'->3'): {hebra_lider}")

    # Hebra rezagada

    print("\n> HEBRA REZAGADA")
    print("  Se forma de manera discontinua mediante fragmentos de Okazaki.")
    print("  Cada fragmento necesita su propio cebador.")

    fragmentos = []

    for i in range(0, len(hebra_5_3), TAM_FRAGMENTO):

        trozo_molde = hebra_5_3[i:i + TAM_FRAGMENTO]
        fragmento_nuevo = complementaria_inversa_adn(trozo_molde)
        k = min(TAM_CEBADOR, len(fragmento_nuevo))
        cebador = cebador_arn(fragmento_nuevo, k)

        fragmentos.append(fragmento_nuevo)

        print(f"\n  Fragmento de Okazaki {len(fragmentos)}:")
        print(f"    Molde: 5'-{trozo_molde}-3'")
        print(f"    PRIMASA -> cebador de ARN: 5'-{cebador}-3'")
        print(f"    ADN POLIMERASA III -> 5'-{cebador}{fragmento_nuevo[k:]}-3'  (ARN + ADN)")
        print(f"    ADN POLIMERASA I -> cebador sustituido: 5'-{fragmento_nuevo}-3'")

    print("\n  La LIGASA une los fragmentos de Okazaki.")

    hebra_rezagada = "".join(reversed(fragmentos))

    print(f"  Hebra rezagada final (5'->3'): {hebra_rezagada}")

    linea()

    print("Resultado de la replicación:")
    print("Se obtienen 2 moleculas de ADN.")
    print("Cada una contiene una hebra original y una hebra nueva.")

    mostrar_doble_helice(hebra_lider, hebra_3_5, "Molecula de ADN 1:")

    print("\nMolecula de ADN 2:")
    print(f"  Hebra original (5'-> 3'): {hebra_5_3}")
    print(f"  Hebra nueva     (5'-> 3'): {hebra_rezagada}")

    return hebra_lider, hebra_rezagada


# 4. TRANSCRIPCIÓN

def transcripcion(hebra_molde):
    """
    Genera el ARNm a partir de la hebra molde de ADN.
    """

    linea("PASO 2: TRANSCRIPCIÓN  (ADN -> ARNm)")

    print(f"Hebra molde (3'-> 5'): {hebra_molde}")

    print("\n> ARN POLIMERASA: utiliza la hebra molde para formar el ARNm.")
    
    arnm = "".join(COMPLEMENTO_ARN[b] for b in hebra_molde)

    print(f"\n  ADN molde (3'-> 5'): {hebra_molde}")
    print(f"  ARNm       (5'-> 3'): {arnm}")

    return arnm


# 5. TRADUCCIÓN

def traduccion(arnm):
    """
    Traduce el ARNm a una secuencia de aminoácidos.
    """

    linea("PASO 3: TRADUCCIÓN  (ARNm -> PROTEINA)")

    print(f"ARNm: {arnm}")

    inicio = arnm.find(CODON_INICIO)

    if inicio == -1:
        print("\nNo se ha encontrado el codon de inicio AUG.")
        return []

    print(f"\n> El ribosoma lee el ARNm en grupos de tres bases llamados codones.")
    print(" Cada codon corresponde a un aminoácido.")

    proteina = []
    i = inicio

    while i + 3 <= len(arnm):

        codon = arnm[i:i + 3]

        if codon in STOP_CODONES:
            print(
                f"  Codon {codon}: Codon de parada. "
                "Finaliza la traduccion."
            )
            break

        aminoacido = CODIGO_GENETICO.get(codon)

        if aminoacido is None:
            print(f"  Codon {codon}: no reconocido.")
            break

        etiqueta_inicio = " (inicio)" if codon == CODON_INICIO else ""

        print(f"  Codon {codon} -> {aminoacido}{etiqueta_inicio}")

        proteina.append(aminoacido)
        i += 3

    else:
        print("  Se ha terminado el ARNm sin encontrar un codon de parada.")

    linea()

    print("Proteina final:")
    print("  " + (" - ".join(proteina) if proteina else "(vacia)"))

    return proteina


# 6. PROGRAMA 


def elegir_origen_adn():

    print("¿Como quieres introducir la molecula de ADN inicial?")
    print("  1: Generarla aleatoriamente")
    print("  2: Escribirla por consola")
    print("  3: Leerla desde un fichero de texto")

    opcion = input("Opcion [1/2/3] (por defecto 1): ").strip()

    if opcion == "2":
        return pedir_cadena_manual()

    if opcion == "3":
        ruta = input("Ruta del fichero: ").strip()
        return leer_cadena_de_fichero(ruta)

    return generar_cadena_aleatoria()


def main():

    linea("SIMULADOR DEL DOGMA CENTRAL DE LA BIOLOGIA MOLECULAR")

    hebra_5_3, hebra_3_5 = elegir_origen_adn()

    hebra_lider, hebra_rezagada = replicacion(
        hebra_5_3,
        hebra_3_5
    )

    print("\nContinuamos con la transcripción de la hebra molde original.")

    arnm = transcripcion(hebra_3_5)

    traduccion(arnm)

    linea("FIN DE LA SIMULACIÓN")


if __name__ == "__main__":
    main()