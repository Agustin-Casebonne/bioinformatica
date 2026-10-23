from Bio import Entrez, SeqIO

Entrez.email = "jaime.rivero103@alu.ulpgc.es"

# Descarga en formato GenBank
handle = Entrez.efetch(db="nucleotide", id="NC_045512",
                       rettype="gb", retmode="text")
genoma = SeqIO.read(handle, "genbank")
handle.close()

print(f"Genoma: {genoma.description}\n")

# Recorrer features, filtrar "gene" y listar nombre y coordenadas
contador = 0
for caracteristica in genoma.features:
    if caracteristica.type == "gene":
        contador += 1
        nombre = caracteristica.qualifiers.get("gene", ["Desconocido"])[0]
        inicio = int(caracteristica.location.start) + 1  # Biopython cuenta desde 0
        fin = int(caracteristica.location.end)
        print(f"{contador}. {nombre}  (inicio: {inicio}, fin: {fin})")

print(f"\nTotal de genes: {contador}")