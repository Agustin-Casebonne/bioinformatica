from Bio import SeqIO


# Leer la secuencia completa descargada de NCBI
registro = SeqIO.read("sequence.fasta", "fasta")
secuencia = registro.seq

# Extraer la CDS indicada por NCBI: posiciones 60-392
cds = secuencia[59:392]

# Guardar la CDS en formato FASTA
with open("adn.fasta", "w") as archivo:
    archivo.write(">INS_CDS_NCBI_NM_000207.3\n")
    archivo.write(str(cds) + "\n")

print("CDS extraída correctamente.")
print("Longitud:", len(cds), "nucleótidos")
print("Secuencia:")
print(cds)