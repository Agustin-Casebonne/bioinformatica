import random

adn_dicc = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A"
}

arn_dicc = {
    "A": "U",
    "C": "G",
    "G": "C",
    "T": "A"
}


codon_a_aminoacido = {
    "UUU": "Fenilalanina", "UUC": "Fenilalanina", "UUA": "Leucina", "UUG": "Leucina",
    "CUU": "Leucina", "CUC": "Leucina", "CUA": "Leucina", "CUG": "Leucina",
    "AUU": "Isoleucina", "AUC": "Isoleucina", "AUA": "Isoleucina", "AUG": "Metionina",
    "GUU": "Valina", "GUC": "Valina", "GUA": "Valina", "GUG": "Valina",
    "UCU": "Serina", "UCC": "Serina", "UCA": "Serina", "UCG": "Serina",
    "CCU": "Prolina", "CCC": "Prolina", "CCA": "Prolina", "CCG": "Prolina",
    "ACU": "Treonina", "ACC": "Treonina", "ACA": "Treonina", "ACG": "Treonina",
    "GCU": "Alanina", "GCC": "Alanina", "GCA": "Alanina", "GCG": "Alanina",
    "UAU": "Tirosina", "UAC": "Tirosina",
    "CAU": "Histidina", "CAC": "Histidina", "CAA": "Glutamina", "CAG": "Glutamina",
    "AAU": "Asparagina", "AAC": "Asparagina", "AAA": "Lisina", "AAG": "Lisina",
    "GAU": "Ácido aspártico", "GAC": "Ácido aspártico", "GAA": "Ácido glutámico", "GAG": "Ácido glutámico",
    "UGU": "Cisteína", "UGC": "Cisteína", "UGG": "Triptófano",
    "CGU": "Arginina", "CGC": "Arginina", "CGA": "Arginina", "CGG": "Arginina",
    "AGU": "Serina", "AGC": "Serina", "AGA": "Arginina", "AGG": "Arginina",
    "GGU": "Glicina", "GGC": "Glicina", "GGA": "Glicina", "GGG": "Glicina",
}

def generar_cadena_adn():
    cadena = ""
    for i in range(9):
        base = random.choice(["A", "C", "G", "T"])
        cadena += base
    cadena_2 = ""
    for i in cadena:
        cadena_2 += adn_dicc[i]
    return cadena, cadena_2

def replicacion(helice1, helice2):
    print("Helicasa: abriendo la doble hélice...")

    helice1_expandida = ""
    print("Hebra líder (síntesis continua, ADN polimerasa, 1 cebador):")
    for i in range(len(helice1)):
        helice1_expandida += adn_dicc[helice1[i]]
    print("  ", helice1, "->", helice1_expandida)

    helice2_expandida = ""
    print("Hebra rezagada (síntesis discontinua, fragmentos de Okazaki):")
    for i in range(0, len(helice2), 3):
        trozo = helice2[i:i+3]
        fragmento = ""
        for base in trozo:
            fragmento += adn_dicc[base]
        helice2_expandida += fragmento
        print(f"  Cebador + fragmento Okazaki: {trozo} -> {fragmento}")
    print("  Ligasa: uniendo fragmentos ->", helice2_expandida)

    return helice1_expandida, helice2_expandida

def transcripcion(helice1, helice2):
    helice1_traducida = ""
    helice2_traducida = ""
    for i in range(len(helice1)):
        helice1_traducida += arn_dicc[helice1[i]]
        helice2_traducida += arn_dicc[helice2[i]]
    print("helice 1 transcrita:", helice1,  helice1_traducida)
    print("helice 2 transcrita:", helice2, helice2_traducida)
    return helice1_traducida, helice2_traducida


def traducir(helice1, helice2):
    helice1_traducida = ""
    helice2_traducida = ""
    codon1 = ""
    codon2 = ""
    stop_codones = ("UAA", "UAG", "UGA")
    parar1 = False
    parar2 = False
    for i in range(0, len(helice1), 3):
        c1 = helice1[i:i+3]
        c2 = helice2[i:i+3]
        if not parar1:
            if c1 in stop_codones:
                print(f"Helice 1: codón {c1} -> señal de terminación")
                parar1 = True
            else:
                helice1_traducida += codon_a_aminoacido[c1] + " "
                codon1 += c1 + " "
        if not parar2:
            if c2 in stop_codones:
                print(f"Helice 2: codón {c2} -> señal de terminación")
                parar2 = True
            else:
                helice2_traducida += codon_a_aminoacido[c2] + " "
                codon2 += c2 + " "
    print("Aminoácidos de la helice 1:", codon1, helice1_traducida)
    print("Aminoácidos de la helice 2:", codon2, helice2_traducida)
    return helice1_traducida, helice2_traducida

if __name__ == "__main__":
    print("Generando cadena de ADN...")
    helice1, helice2 = generar_cadena_adn()
    print("Cadena de ADN generada:", helice1, helice2)
    print("Replicando: Generando divición de helices...")
    helice1_replicada, helice2_replicada = replicacion(helice1, helice2)
    print("Transcribiendo: Generando ARN...")
    helice1_transcrita, helice2_transcrita = transcripcion(helice1_replicada, helice2_replicada)
    traducir(helice1_transcrita, helice2_transcrita)