from Bio.Seq import Seq


# Secuencia de ARNm
arnm = Seq("AUGUAUGCUUAA")

# Traducción del ARNm a proteína
proteina = arnm.translate()

print("ARNm:")
print("5'-" + str(arnm) + "-3'")

print("\nProteína:")
print(proteina)

print("\nResultado esperado:")
print("Metionina - Tirosina - Alanina")