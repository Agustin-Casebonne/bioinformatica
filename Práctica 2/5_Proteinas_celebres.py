from Bio import Entrez, SeqIO

Entrez.email = "jaime.rivero103@alu.ulpgc.es"

proteinas = {
    "Hemoglobina (subunidad alfa)": ("P69905",
        "Transporta oxígeno desde los pulmones a los tejidos en los glóbulos rojos."),
    "Mioglobina": ("P02144",
        "Almacena oxígeno en el músculo y lo libera durante el esfuerzo."),
    "Oxitocina (precursor)": ("P01178",
        "Hormona que participa en las contracciones del parto, la lactancia y el vínculo social."),
}

informe = []

for nombre, (accesion, funcion) in proteinas.items():
    print("=" * 70)
    print(f"{nombre}  ->  {accesion}")

    # Información en formato GenPept
    handle = Entrez.efetch(db="protein", id=accesion, rettype="gp", retmode="text")
    gp = SeqIO.read(handle, "genbank")
    handle.close()

    organismo = gp.annotations.get("organism", "?")
    print("Descripción:", gp.description)
    print("Organismo:  ", organismo)
    print("Longitud:   ", len(gp.seq), "aminoácidos")
    print("Función:    ", funcion)

    # Descarga de archivo FASTA y guardado en disco
    handle = Entrez.efetch(db="protein", id=accesion, rettype="fasta", retmode="text")
    fasta = SeqIO.read(handle, "fasta")
    handle.close()

    fichero = f"{accesion}.fasta"
    SeqIO.write(fasta, fichero, "fasta")
    print(f"FASTA guardado en: {fichero}")

    informe.append(
        f"{nombre} ({accesion})\n"
        f"  Descripción: {gp.description}\n"
        f"  Organismo: {organismo}\n"
        f"  Longitud: {len(gp.seq)} aminoácidos\n"
        f"  Función: {funcion}\n"
        f"  Archivo FASTA: {fichero}\n"
    )

with open("informe_proteinas.txt", "w", encoding="utf-8") as f:
    f.write("\n".join(informe))
print("\nInforme guardado en 'informe_proteinas.txt'.")