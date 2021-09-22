'''
NAME
    ProteinDB.py
VERSION
    [1.0]
AUTHOR
    Daianna Gonzalez Padilla <daianna@lcg.unam.mx>
DESCRIPTION
    This programs takes a PDB file and returns certain residues of a protein chain 
CATEGORY
     Protein Data Bank files analysis
USAGE
    None
ARGUMENTS
    This program doesn't take arguments
INPUT
    The input of the function is the path to the file, the chain name and the residue name
OUTPUT
    Returns a list whose elements are lists with the residue name and its ID
    
EXAMPLES
    Example 1: gets  get_them=get_residue('C:/Users/hp/Downloads/1kcw.pdb','A', 'CYS')
                     print(get_them)
               and returns  [['CYS', 155], ['CYS', 181], ['CYS', 221], ['CYS', 257], ['CYS', 319]
                             ['CYS', 338], ['CYS', 515], ['CYS', 541], ['CYS', 618], ['CYS', 680],
                             ['CYS', 699], ['CYS', 855], ['CYS', 881], ['CYS', 1021]]
SOURCE
    https://github.com/daianna21/python_class/blob/master/tareas/ProteinDB.py
'''



from Bio import PDB

def get_residue(path, chain_name, res_name):
    """
    This function gets a PDB file, the name of a protein's module chain and a certain residue
    and returns the IDs of those residues in the specified chain
            Parameters:
                    path (str): absolute path to file
                    chain_name (str): Single letter name of the protein module chain
                    res_name (str): Three leter name of the residue to search
            Returns:
                    residues (list): list of lists with the residue name and the ID of each one
    """
    # List to save the residues
    residues=[]
    # Create the parser and ignore warnings
    parser = PDB.PDBParser(QUIET=True)
    # Create the structure from the file 
    struct=parser.get_structure('protein',path)
    # For each module, get each chain and search the given one 
    for model in struct:
        for chain in model:
            if chain==model[chain_name]:
                # Search the residues in the chain and store them with its IDs
                for residue in chain:
                    if residue.get_resname() == res_name:
                        residues.append([residue.get_resname(), residue.get_id()[1]])
                        #print(residue)
    return(residues)

    

#Example
get_them=get_residue('C:/Users/hp/Downloads/1kcw.pdb','A', 'THR')
print(get_them)
 

