"""
Ejercicio 4: A la caza de proteínas de resistencia
Plásmido pT181 de Staphylococcus aureus (NC_001371).
Descarga el GenBank, extrae todas las CDS y las traduce a proteína.
"""
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord

Entrez.email = "tu_correo@ejemplo.com"
ACCESO = "NC_001371"

# "gbwithparts" en vez de "gb": con "gb" el NCBI a veces devuelve el registro sin secuencia
with Entrez.efetch(db="nucleotide", id=ACCESO, rettype="gbwithparts", retmode="text") as handle:
    plasmido = SeqIO.read(handle, "genbank")

print(f"{plasmido.id} - {plasmido.description}")
print(f"Longitud: {len(plasmido.seq)} pb\n")

proteinas = []
for feature in plasmido.features:
    if feature.type != "CDS":
        continue

    q = feature.qualifiers
    gen = q.get("gene", q.get("locus_tag", ["?"]))[0]
    producto = q.get("product", ["producto desconocido"])[0]
    tabla = int(q.get("transl_table", [11])[0])   # 11 = código genético bacteriano

    # extract() respeta la hebra (complemento inverso si strand = -1)
    # y las CDS que cruzan el origen del plásmido circular
    adn_cds = feature.extract(plasmido.seq)
    proteina = adn_cds.translate(table=tabla, to_stop=True)

    # Comprobación con la traducción oficial del NCBI, si existe
    oficial = q.get("translation", [None])[0]
    coincide = "sí" if oficial == str(proteina) else ("no" if oficial else "n/d")

    print(f"Gen: {gen} | Producto: {producto}")
    print(f"  Posición: {feature.location}")
    print(f"  Longitud: {len(proteina)} aa | Coincide con NCBI: {coincide}")
    print(f"  {proteina}\n")

    proteinas.append(SeqRecord(proteina, id=f"{gen}", description=f"{producto} [{ACCESO}]"))

print(f"Total de CDS traducidas: {len(proteinas)}")
SeqIO.write(proteinas, "proteinas_pT181.fasta", "fasta")
print("Proteínas guardadas en 'proteinas_pT181.fasta'")
