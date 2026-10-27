from Bio import Entrez, SeqIO

Entrez.email = "tu_correo@ejemplo.com"   # el NCBI exige un correo identificativo

TERMINO = "Measles morbillivirus[Organism] AND refseq[filter]"

# Búsqueda en la base de datos 'nucleotide'
with Entrez.esearch(db="nucleotide", term=TERMINO) as handle:
    resultados = Entrez.read(handle)

print(f"Resultados encontrados: {resultados['Count']}")
if not resultados["IdList"]:
    raise SystemExit("La búsqueda no devolvió ningún ID.")

# Primer ID
id_genoma = resultados["IdList"][0]
print(f"Primer ID: {id_genoma}")

# Descarga en formato FASTA
with Entrez.efetch(db="nucleotide", id=id_genoma, rettype="fasta", retmode="text") as handle:
    registro = SeqIO.read(handle, "fasta")

# 4) Descripción y primeros 100 nucleótidos
print(f"\nDescripción: {registro.description}")
print(f"Longitud total: {len(registro.seq)} nt")
print(f"Primeros 100 nucleótidos:\n{registro.seq[:100]}")

# Opcional: guardar el genoma
SeqIO.write(registro, "sarampion.fasta", "fasta")
