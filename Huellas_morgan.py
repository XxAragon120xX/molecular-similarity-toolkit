import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import AllChem, Draw

# Función para calcular y mostrar las huellas Morgan
def mostrar_huellas_morgan(smiles):
    # Convertir la cadena SMILES a una molécula RDKit
    molecula = Chem.MolFromSmiles(smiles)
    
    # Verificar si la molécula se generó correctamente
    if not molecula:
        print("SMILES no válida. Inténtalo de nuevo.")
        return
    
    # Calcular las huellas Morgan
    onbits = {}
    huella_morgan = AllChem.GetMorganFingerprintAsBitVect(molecula, 2, nBits=512, bitInfo=onbits)
    
    # Mostrar los índices de los bits activos
    print("Índices de bits activos:", tuple(huella_morgan.GetOnBits()))
    print("Número de bits activos:", len(huella_morgan.GetOnBits()))
    
    # Dibujar los fragmentos responsables de los bits activos
    #for bit in huella_morgan.GetOnBits():
    #    img = Draw.DrawMorganBit(molecula, bit, onbits)
    #    plt.imshow(img)
    #    plt.title(f"Fragmento responsable del bit {bit}")
    #    plt.axis('off')
    #    plt.show()

# Ingresar la cadena SMILES
smiles = input("Introduce una cadena SMILES: ")
mostrar_huellas_morgan(smiles)

#//////////////////////////////////////////////////////////////////////////////

def calcular_similitud(lista1, lista2):
    # Convertir las listas en conjuntos para eliminar duplicados
    set1 = set(lista1)
    set2 = set(lista2)

    # Encontrar la intersección entre ambos conjuntos
    interseccion = set1.intersection(set2)

    # Calcular el porcentaje de similitud
    if len(set1) + len(set2) == 0:
        return 0  # Si ambas listas están vacías
    porcentaje_similitud = (len(interseccion) / ((len(set1) + len(set2)) / 2)) * 100

    return porcentaje_similitud

# Solicitar al usuario ingresar las listas
lista1 = input("Introduce la primera lista de números separados por comas: ")
lista2 = input("Introduce la segunda lista de números separados por comas: ")

# Convertir las cadenas de entrada a listas de números
lista1 = [int(x) for x in lista1.split(",")]
lista2 = [int(x) for x in lista2.split(",")]

# Calcular la similitud
similitud = calcular_similitud(lista1, lista2)

# Mostrar el resultado
print(f"La similitud entre las dos listas es de {similitud:.2f}%")
