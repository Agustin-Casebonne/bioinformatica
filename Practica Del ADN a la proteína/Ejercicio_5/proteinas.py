from Bio.Seq import Seq


# Secuencia de aminoácidos
proteina = Seq("MISGVKH")


print("Secuencia de la proteína:")
print(proteina)

print("\nExtremo N-terminal:")
print(proteina[0])

print("\nExtremo C-terminal:")
print(proteina[-1])

print("\nLongitud:")
print(len(proteina), "aminoácidos")