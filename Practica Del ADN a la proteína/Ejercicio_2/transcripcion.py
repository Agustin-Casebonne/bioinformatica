from Bio.Seq import Seq
from Bio import SeqIO


# Leer la secuencia de ADN desde un archivo FASTA
archivo_fasta = SeqIO.read("adn.fasta", "fasta")
adn_molde = archivo_fasta.seq

# Transcripción del ADN molde a ARN
arnm = adn_molde.complement().transcribe()

print("ADN molde:")
print("3'-" + str(adn_molde) + "-5'")

print("\nARNm:")
print("5'-" + str(arnm) + "-3'")

print("\nResultado esperado:")
print("5'-AUGCCUGAAUGC-3'")

print("\n¿Coinciden los resultados?")

if str(arnm) == "AUGCCUGAAUGC":
    print("Sí, el resultado coincide.")
else:
    print("No, los resultados no coinciden.")