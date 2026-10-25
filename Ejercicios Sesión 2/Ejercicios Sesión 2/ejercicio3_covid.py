"""
Ejercicio 3: La COVID-19
Lista los genes del genoma de referencia del SARS-CoV-2 (NC_045512).
"""
import os
from Bio import Entrez, SeqIO
from dotenv import load_dotenv

# El NCBI exige un correo identificativo; se lee del fichero .env (no se sube a GitHub)
load_dotenv()
Entrez.email = os.getenv("NCBI_EMAIL")
if not Entrez.email:
    raise SystemExit("Falta NCBI_EMAIL: añádelo al fichero .env de la raíz del proyecto")
ACCESO = "NC_045512"

# Descarga en formato GenBank
with Entrez.efetch(db="nucleotide", id=ACCESO, rettype="gb", retmode="text") as handle:
    registro = SeqIO.read(handle, "genbank")

print(f"{registro.id} - {registro.description}\n")

# Recorrer features, quedarse solo con los de tipo 'gene' y sacar su nombre
genes = []
for feature in registro.features:
    if feature.type == "gene":
        nombre = feature.qualifiers.get("gene", ["(sin nombre)"])[0]
        inicio = int(feature.location.start) + 1   # Biopython usa coordenadas 0-based
        fin = int(feature.location.end)
        hebra = "+" if feature.location.strand == 1 else "-"
        genes.append((nombre, inicio, fin, hebra))

# Listado numerado
print(f"Número de genes: {len(genes)}\n")
print(f"{'Nº':<4}{'Gen':<10}{'Inicio':>8}{'Fin':>8}{'Long.':>8}  Hebra")
for i, (nombre, inicio, fin, hebra) in enumerate(genes, start=1):
    print(f"{i:<4}{nombre:<10}{inicio:>8}{fin:>8}{fin - inicio + 1:>8}  {hebra}")
