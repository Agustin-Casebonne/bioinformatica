# Ejercicio 5: Tres proteínas "célebres"

Las secuencias proceden de UniProtKB/Swiss-Prot (entradas revisadas, *Homo sapiens*). El script `ejercicio5_proteinas.py` descarga los FASTA desde `https://rest.uniprot.org/uniprotkb/<acceso>.fasta`.

## 1. Hemoglobina (cadena beta) – UniProt P68871 (HBB_HUMAN)
- **Gen:** *HBB*, cromosoma 11.
- **Tamaño:** 147 aa en el precursor (146 en la proteína madura, porque se elimina la Met inicial).
- **Función:** transporta O₂ de los pulmones a los tejidos dentro de los eritrocitos. La hemoglobina adulta (HbA) es un tetrámero α₂β₂, y cada subunidad lleva un grupo hemo con un ion Fe²⁺.
- **Interés:** la mutación Glu→Val en la posición 6 de la cadena beta (HbS) causa la anemia falciforme. Las mutaciones que reducen su síntesis causan la β-talasemia.

## 2. Mioglobina – UniProt P02144 (MYG_HUMAN)
- **Gen:** *MB*, cromosoma 22.
- **Tamaño:** 154 aa.
- **Función:** almacena O₂ en el músculo esquelético y cardíaco y lo cede cuando la demanda es alta. Es un monómero con un solo grupo hemo, y su afinidad por el O₂ es mayor que la de la hemoglobina.
- **Interés:** fue la primera proteína cuya estructura 3D se resolvió por cristalografía de rayos X (Kendrew, 1958; Nobel 1962, con Perutz por la hemoglobina). En clínica se usa como marcador de daño muscular.

## 3. Oxitocina – UniProt P01178 (NEU1_HUMAN)
- **Gen:** *OXT*, cromosoma 20.
- **Tamaño:** el precursor (oxitocina-neurofisina 1) tiene 125 aa. La hormona madura es un nonapéptido, CYIQNCPLG-NH₂, con un puente disulfuro entre las dos cisteínas.
- **Función:** se sintetiza en el hipotálamo y se libera desde la neurohipófisis. Provoca las contracciones uterinas en el parto y la eyección de leche en la lactancia, e interviene en el vínculo social y afectivo.
- **Interés:** demuestra que una hormona pequeña se produce cortando un precursor mayor. Su síntesis química (du Vigneaud) le valió el Nobel de Química de 1955.

## Formato FASTA (ejemplo de cabecera)
```
>sp|P68871|HBB_HUMAN Hemoglobin subunit beta OS=Homo sapiens OX=9606 GN=HBB PE=1 SV=2
```
`sp` indica Swiss-Prot. Le siguen el código de acceso y el nombre de la entrada, y después el organismo (OS), el taxón (OX), el gen (GN), el nivel de evidencia (PE) y la versión de la secuencia (SV).
