"""
Ejercicio 5: Proteínas "célebres"
Descarga desde UniProt (Swiss-Prot, revisadas) los FASTA de tres proteínas
humanas de la lista: hemoglobina (cadena beta), mioglobina y oxitocina.
"""
import io
import urllib.request
from Bio import SeqIO

PROTEINAS = {
    "hemoglobina_beta": "P68871",   # HBB_HUMAN
    "mioglobina": "P02144",         # MYG_HUMAN
    "oxitocina": "P01178",          # NEU1_HUMAN (precursor oxitocina-neurofisina 1)
}

registros = []
for nombre, acceso in PROTEINAS.items():
    url = f"https://rest.uniprot.org/uniprotkb/{acceso}.fasta"
    with urllib.request.urlopen(url) as resp:
        texto = resp.read().decode()

    with open(f"{nombre}.fasta", "w") as f:
        f.write(texto)

    rec = SeqIO.read(io.StringIO(texto), "fasta")
    registros.append(rec)
    print(f"{nombre:<18} {rec.id:<25} {len(rec.seq):>4} aa")
    print(f"  {rec.description}\n")

SeqIO.write(registros, "tres_proteinas.fasta", "fasta")
print("FASTA combinado: tres_proteinas.fasta")
