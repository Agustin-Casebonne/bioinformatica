from Bio import Entrez, SeqIO

Entrez.email = "jaime.rivero103@alu.ulpgc.es"

# Búsqueda en la base secundaria (RefSeq)
consulta = "Measles morbillivirus[Organism] AND srcdb_refseq[Properties]"
handle = Entrez.esearch(db="nucleotide", term=consulta)
resultado = Entrez.read(handle)
handle.close()

# Primer ID
if resultado["IdList"]:
    id_objetivo = resultado["IdList"][0]
    print("ID encontrado:", id_objetivo)

    # Descarga FASTA
    handle = Entrez.efetch(db="nucleotide", id=id_objetivo,
                           rettype="fasta", retmode="text")
    registro = SeqIO.read(handle, "fasta")
    handle.close()

    # Descripción y primeros 100 nucleótidos
    print("Descripción:", registro.description)
    print("Primeros 100 nucleótidos:", registro.seq[:100])
else:
    print("No se encontraron resultados.")