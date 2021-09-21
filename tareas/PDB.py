'''
NAME
    PDB.py
VERSION
    [1.0]
AUTHOR
    Daianna Gonzalez Padilla <daianna@lcg.unam.mx>
DESCRIPTION
    This programs takes a PDB file and returns which Cys pairs of a protein chain may form disulfide bonds
CATEGORY
     Protein Data Bank files analysis
USAGE
    None
ARGUMENTS
    This program doesn't take arguments
INPUT
    The input of the function is the path to the file and the chain name
OUTPUT
    Returns the cisteins ID's and the distance between their sulfur atoms
    
EXAMPLES
    Example 1: gets  bonds=ds_bond('C:/Users/hp/Downloads/1kcw.pdb','A')
               and returns  Possible disulfide bond: Cys  155  --> Cys  181  Distance:  1.9980831
                            Possible disulfide bond: Cys  257  --> Cys  338  Distance:  2.0317485
                            Possible disulfide bond: Cys  515  --> Cys  541  Distance:  2.0152018
                            Possible disulfide bond: Cys  618  --> Cys  699  Distance:  2.027114
                            Possible disulfide bond: Cys  855  --> Cys  881  Distance:  2.0152495
SOURCE
    https://github.com/daianna21/python_class/blob/master/tareas/PDB.py
'''



from Bio import PDB

def ds_bond(path, chain_name):
    """
    This function gets a PDB file and the name of a protein's module chain  and returns
    possible disulfide bonds within the chain
            Parameters:
                    path (str): absolute path to file
                    chain_name (str): Single letter name of the protein module chain 
            Returns:
                    ID of cisteines involved in each bond and the distance between their
                    sulfur atoms in Armstrongs
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

 

