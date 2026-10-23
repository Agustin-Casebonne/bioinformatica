"""
Ejercicio 3: La COVID-19
Lista los genes del genoma de referencia del SARS-CoV-2 (NC_045512).
"""
from Bio import Entrez, SeqIO

Entrez.email = ""
ACCESO = "NC_045512"

# 1) Descarga en formato GenBank
with Entrez.efetch(db="nucleotide", id=ACCESO, rettype="gb", retmode="text") as handle:
    registro = SeqIO.read(handle, "genbank")

print(f"{registro.id} - {registro.description}\n")

# 2-4) Recorrer features, quedarse solo con los de tipo 'gene' y sacar su nombre
genes = []
for feature in registro.features:
    if feature.type == "gene":
        nombre = feature.qualifiers.get("gene", ["(sin nombre)"])[0]
        inicio = int(feature.location.start) + 1   # Biopython usa coordenadas 0-based
        fin = int(feature.location.end)
        hebra = "+" if feature.location.strand == 1 else "-"
        genes.append((nombre, inicio, fin, hebra))

# 5) Listado numerado
print(f"Número de genes: {len(genes)}\n")
print(f"{'Nº':<4}{'Gen':<10}{'Inicio':>8}{'Fin':>8}{'Long.':>8}  Hebra")
for i, (nombre, inicio, fin, hebra) in enumerate(genes, start=1):
    print(f"{i:<4}{nombre:<10}{inicio:>8}{fin:>8}{fin - inicio + 1:>8}  {hebra}")
