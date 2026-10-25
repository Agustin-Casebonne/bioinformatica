from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

proteinas = []

# Recorremos cada secuencia del fichero FASTA
for registro in SeqIO.parse("secuencias.fasta", "fasta"):

    # Porcentaje de GC
    gc = gc_fraction(registro.seq) * 100
    print(registro.id, "-> GC:", round(gc, 2), "%")

    # Traducir ADN a aminoácidos
    registro.seq = registro.seq.translate()
    print("   Proteína:", registro.seq)

    proteinas.append(registro)

# Guardar las proteínas en un FASTA nuevo (mismo id y descripción)
SeqIO.write(proteinas, "proteinas.fasta", "fasta")
