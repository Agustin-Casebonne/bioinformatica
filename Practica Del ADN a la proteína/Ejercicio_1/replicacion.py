from Bio.Seq import Seq


# Secuencia de ADN original
adn = Seq("ATGCCGTTAGCT")

# Generar la hebra complementaria
complementaria = adn.complement()


print("ADN original:")
print("5'-" + str(adn) + "-3'")

print("\nHebra complementaria:")
print("3'-" + str(complementaria) + "-5'")

print("\nResultado manual:")
print("3'-TACGGCAATCGA-5'")

print("\n¿Coinciden los resultados?")

if str(complementaria) == "TACGGCAATCGA":
    print("Sí, el resultado coincide.")
else:
    print("No, los resultados no coinciden.")