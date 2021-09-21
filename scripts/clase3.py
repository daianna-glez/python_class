# Clase 3: Estructura 3D de las proteínas
## Formato PDB
'''
Base de datos con las estructuras 3D de las proteínas y ácidos
nucleicos

Se crean objetos structure tipo SMCRA

'''
# ESTRUCTURA S
#Para crear un structure
from Bio import PDB 
#QUIET=True para no ver los warnings, crear parser
parser = PDB.PDBParser(QUIET=True)

struc=parser.get_structure('prot','C:/Users/hp/Downloads/1fat.pdb')
struc.child_dict

#Ver cuántos y cúales modelos tiene el objeto
print(struc.child_dict)
print(struc.child_list)

#Acceder a metadatos con el atributo header

print(struc.header.keys())
print(struc.header['structure_method'])
print(struc.header['resolution'])
print(struc.header['source'])
print(struc.header['journal'])


#Ejericio 1:
parser = PDB.PDBParser(QUIET=True)
struct=parser.get_structure('prot','C:/Users/hp/Downloads/1kcw.pdb')
print(struct.header['structure_method'])
print(struct.header['resolution'])


#MODELOS M
## Para acceder a los modelos podemos iterar en el objeto structure

modelo = struct[0]
print(modelo)
for model in struc:
    print(model)
    print(model.child_dict)


#CADENAS C

for chain in model:
    print(chain)
    #print(chain.child_dict)
## O con chain = model['A'] o chain = model.child_list[0]


#RESIDUO R

#for residue in chain:
  #print(residue)
residue = chain[1]
print(residue)
## Obtener ID del residuo (segundo elemento de get_id):
print(residue.get_id())
print(residue.get_id()[1])
  
##Obtener nombre del residuo en 3 letras:
print(residue.get_resname())

#Ubicar un residuo que nos interesa:
for residue in chain:
  if residue.get_resname() == "PHE":
    print(residue)

#Ubicar los àtomos con child_list o dict
residue = chain[1]
print(residue.child_list)




 #Ejercicio 2:
cys_A=[]
for model in struct:
    for chain in model:
        if chain==model['A']:
            for residue in chain:
                if residue.get_resname() == 'CYS':
                    cys_A.append(residue)

print(cys_A)





#ATOMO
residue = chain[22]
for atom in residue:
    print(atom)# o con atom = residue['CA'] o atom = residue.child_list[1]
##Información de su posición 3D, elemento que es y ID
    #print(atom.coord)
    #print(atom.element)
    #print(atom.id)

## Distancia entre dos átomos:
atomo_1 = residue['CA']
atomo_2 = residue['N']
print(atomo_1-atomo_2)




#Ejercicio 3:
print(cys_A[0])
#Cada atomo de un residuo 
for atom in cys_A[0]:
    print(atom.element, atom.id)

#Ejercicio 4:
pares=[]    
for cys_1 in cys_A:
    for cys_2 in cys_A:
        if (cys_1['SG']-cys_2['SG'])<8:
            pares.append(cys_1['SG']-cys_2['SG'])
print(pares)


#Tarea:
pares=[]    
for cys_1 in cys_A:
        for cys_2 in cys_A:
            if not cys_1==cys_2:
                if (cys_1['SG']-cys_2['SG'])<8 and (cys_1['SG']-cys_2['SG']) not in pares:
                    pares.append(cys_1['SG']-cys_2['SG'])
                    print(cys_1.get_id()[1], cys_2.get_id()[1],cys_1['SG']-cys_2['SG'])
#print(pares)
                    
from Bio import PDB

def ds_bond(path, chain_name):
    """
    This function gets a PDB file and the name of a protein's module chain  and returns
    possible disulfide bonds within the chain
            Parameters:
                    path (str): absolute path to file
                    chain_name (str): Single letter name of the protein's module chain 
            Returns:
                    ID of cisteines involved in each bond and the distance between them
                    in Armstrongs
    """
    # List to save Cys residues
    pairs=[]
    cys_A=[]
    # Create the parser and ignore warnings
    parser = PDB.PDBParser(QUIET=True)
    # Create the structure from the file 
    struct=parser.get_structure('protein',path)
    # For each module, get each chain and search the given one 
    for model in struct:
        for chain in model:
            if chain==model[chain_name]:
                # Search Cys residues in the chain and store them 
                for residue in chain:
                    if residue.get_resname() == 'CYS':
                        cys_A.append(residue)
    # For each unique pair of Cys residues get its ID's and its sulfur atoms                     
    for cys_1 in cys_A:
            for cys_2 in cys_A:
                if not cys_1==cys_2:
                    # Distance between sulfur atoms less than 8 A for pairs of Cys not repeated
                    if (cys_1['SG']-cys_2['SG'])<8 and [cys_1.get_id()[1], cys_2.get_id()[1]] not in pairs and [cys_2.get_id()[1], cys_1.get_id()[1]] not in pairs:
                        pairs.append([cys_1.get_id()[1], cys_2.get_id()[1]])
                        print("Possible disulfide bond: Cys ",cys_1.get_id()[1], " --> Cys ", cys_2.get_id()[1],
                              " Distance: ",cys_1['SG']-cys_2['SG'])

    

#Example
bonds=ds_bond('C:/Users/hp/Downloads/1kcw.pdb','A')

 

