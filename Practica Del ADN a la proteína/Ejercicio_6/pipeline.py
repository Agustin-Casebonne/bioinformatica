from Bio import SeqIO


# Leer la secuencia de ADN desde un archivo FASTA

archivo_fasta = SeqIO.read("adn.fasta", "fasta")
adn = archivo_fasta.seq


# 1. Replicación
complementaria = adn.complement()


print("ADN original:")
print("5'-" + str(adn) + "-3'")

print("\nHebra complementaria:")
print("3'-" + str(complementaria) + "-5'")


print("\nDespués de la replicación:")

print("\nMolécula 1:")
print("5'-" + str(adn) + "-3' (original)")
print("3'-" + str(complementaria) + "-5' (nueva)")

print("\nMolécula 2:")
print("5'-" + str(adn) + "-3' (nueva)")
print("3'-" + str(complementaria) + "-5' (original)")


# 2. Transcripción
# La hebra complementaria está orientada 3' -> 5',
# por lo que actúa como hebra molde.
adn_molde = complementaria

arnm = adn_molde.complement().transcribe()


print("\nARNm:")
print("5'-" + str(arnm) + "-3'")


# 3. Traducción
proteina = arnm.translate()


print("\nProteína:")
print(proteina)


# 4. Información del proceso
print("\nResumen del proceso:")
print("ADN -> replicación -> ADN")
print("ADN molde -> transcripción -> ARNm")
print("ARNm -> traducción -> proteína")