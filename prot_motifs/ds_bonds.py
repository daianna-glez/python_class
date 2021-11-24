from Bio import PDB
import argparse
import os
import re

# Create the parser
my_parser = argparse.ArgumentParser(description="Program that obtains protein motifs")
# Add an argument to request the input file
my_parser.add_argument("-i", "--input", nargs='+',
                    type=list,
                    help="List of the proteins files paths",
                    required=True)
# Execute the parse_args() method
args = my_parser.parse_args()

def ds_bond(paths):
    
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
    # Get disulfide bonds in each protein/file
    for path in paths:
        print(path)
        # Parse protein name from the path file
        prot_name=path.split('/')[2]
        prot_name=str(prot_name).split('.')[0]
        print(prot_name)

        # Create output dictionaty
        bonds=[]
        c=0
        
        # Create the PDB parser and ignore warnings
        parser = PDB.PDBParser(QUIET=True)
        # Get the structure from the file 
        struct=parser.get_structure(f'protein{prot_name}',path)
        # For each module, get each chain and search the given one 
        for model in struct:
            for chain in model:
                cys_chain=[]
                pairs=[]
    
                # Search Cys residues in the chain and store them 
                for residue in chain:
                    if residue.get_resname() == 'CYS':
                        cys_chain.append(residue)

                # For each unique pair of Cys residues get its ID's and its sulfur atoms                     
                for cys_1 in cys_chain:
                    for cys_2 in cys_chain:
                        if not cys_1==cys_2:
                            # Distance between sulfur atoms less than 8 A for pairs of Cys not repeated
                            if (cys_1['SG']-cys_2['SG'])<8 and [cys_1.get_id()[1], cys_2.get_id()[1]] not in pairs and [cys_2.get_id()[1], cys_1.get_id()[1]] not in pairs:
                                pairs.append([cys_1.get_id()[1], cys_2.get_id()[1]])
                                bonds.append([cys_1.get_id()[1], cys_2.get_id()[1], cys_1['SG']-cys_2['SG'], model, chain])
                                c=c+1
                                
        #print(bonds)
        prot_dic={'name':prot_name, 'num_bonds': c, 'di_bonds':bonds}
        print(prot_dic)
    

#Example
#bondss=ds_bond(['../docs/1kcw.pdb', '../docs/1fat.pdb'])


#Call the function with the arguments given
ds_bond(args.input)
