from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction

fichero_entrada = "insulina_humana_multiple.fasta"
fichero_salida = "proteinas.fasta"

proteinas = []

for registro in SeqIO.parse(fichero_entrada, "fasta"):
    
    # Contenido de GC en porcentaje, mostrado por pantalla
    porcentaje_gc = gc_fraction(registro.seq) * 100
    print(f"{registro.id}: GC = {porcentaje_gc:.2f} %")

    # Traducción a aminoácidos (se recorta a múltiplo de 3 para evitar avisos)
    resto = len(registro.seq) % 3
    secuencia = registro.seq[:len(registro.seq) - resto]
    secuencia_proteica = secuencia.translate()

    # Registro con el identificador y la descripción originales
    proteinas.append(SeqRecord(secuencia_proteica,
                               id=registro.id,
                               description=registro.description))

# Escritura del archivo FASTA con las proteínas
n = SeqIO.write(proteinas, fichero_salida, "fasta")
print(f"\nSe han escrito {n} secuencias de aminoácidos en '{fichero_salida}'.")