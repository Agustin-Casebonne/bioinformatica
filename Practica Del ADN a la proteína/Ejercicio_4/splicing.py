from Bio.Seq import Seq


# Exones del gen
exon1 = Seq("ATGAAA")
exon2 = Seq("CCCGGG")
exon3 = Seq("TTTCCC")
exon4 = Seq("AAATTT")
exon5 = Seq("GGG")


# Splicing alternativo
transcrito_1 = exon1 + exon2 + exon4 + exon5
transcrito_2 = exon1 + exon3 + exon5


print("Transcrito 1 (exones 1-2-4-5):")
print(transcrito_1)

print("\nTranscrito 2 (exones 1-3-5):")
print(transcrito_2)

print("\nLongitud de los transcritos:")
print("Transcrito 1:", len(transcrito_1), "nucleótidos")
print("Transcrito 2:", len(transcrito_2), "nucleótidos")


# Traducción para observar las diferencias en las proteínas
proteina_1 = transcrito_1.translate()
proteina_2 = transcrito_2.translate()

print("\nProteína del transcrito 1:")
print(proteina_1)

print("\nProteína del transcrito 2:")
print(proteina_2)