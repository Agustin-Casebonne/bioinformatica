from Bio import Entrez, SeqIO

Entrez.email = "jaime.rivero103@alu.ulpgc.es"

# El enunciado indica NC_001371, pero ese código es el plásmido ColE1 de E. coli.
# El plásmido pT181 de Staphylococcus aureus es NC_001393.
accesion = "NC_001393"

print(f"Descargando el plásmido {accesion} desde el NCBI...\n")
handle = Entrez.efetch(db="nucleotide", id=accesion,
                       rettype="gbwithparts", retmode="text")
plasmido = SeqIO.read(handle, "genbank")
handle.close()

# Si el registro viniera sin secuencia, la descargamos en FASTA
try:
    str(plasmido.seq[:1])
except Exception:
    print("El registro GenBank no incluía la secuencia; descargando FASTA...")
    handle = Entrez.efetch(db="nucleotide", id=accesion,
                           rettype="fasta", retmode="text")
    plasmido.seq = SeqIO.read(handle, "fasta").seq
    handle.close()

print(f"Analizando: {plasmido.description}")
print("-" * 60)

contador = 0
for caracteristica in plasmido.features:
    if caracteristica.type == "CDS":
        contador += 1
        nombre = caracteristica.qualifiers.get("product", ["Desconocida"])[0]

        # Extraer el gen (respeta la hebra) y traducir con el código bacteriano
        gen = caracteristica.extract(plasmido.seq)
        proteina = gen.translate(table="Bacterial", to_stop=True)

        print(f"Proteína {contador}: {nombre}")
        print(f"Longitud:  {len(proteina)} aminoácidos")
        print(f"Secuencia: {proteina}\n")

print(f"Total de CDS traducidas: {contador}")