'''
NAME
    motifs.py
VERSION
    [2.1.0]
AUTHOR
    Daianna Gonzalez Padilla <daianna@lcg.unam.mx>
DESCRIPTION
    This module is part of prot_motifs package ..
CATEGORY
     PDB files analysis
USAGE
    None
ARGUMENTS

INPUT
    The input depends on the fuction used. Both requiere the database name to search information

    must be in a folder outside the current directory
OUTPUT
    The output depends on the function

EXAMPLES
    Example 1.1: gets  dictionary={"FieldList":'ECNO', "LinkList":"protein_protein_small_genome"}
                       get_descriptions('protein',dictionary)
               and returns:  EC number for enzyme or CAS registry number
                             All proteins from this genome

    Example 1.2: gets  search_articles(data_base="pubmed") whit the following inputs:
                        name= Valeria Mateo-Estrada
                        word_title= Genomic
                        word_title= DNA
                and returns an output file with the following text:
                        Results:
                        1 article(s) found, whose ID(s) (for the first 20) are :['34282943']


SOURCE
    https://github.com/daianna21/python_class/blob/master/tareas/Entrez_search.py
'''
from Bio import PDB
import argparse
import re

# Create the argument parser
my_parser = argparse.ArgumentParser(description="Program that obtains protein motifs")
# Add arguments: choose which motifs to search
# Request input files
my_parser.add_argument("-i", "--input",  nargs='+',
                    type=str,
                    help="List of the proteins files paths",
                    required=True)
# Give distance between S atoms
my_parser.add_argument("-d", "--disulfide",
                    type=float,
                    help="Distance between S-S atoms. 8 by default",
                    required=False)
# Give a pattern to search
my_parser.add_argument("-a", "--alpha",
                    type=str,
                    help="Sequence pattern to search alpha helix. Use default..",
                    required=False)
my_parser.add_argument("-b", "--beta",
                    type=str,
                    help="Sequence pattern to search beta sheets",
                    required=False)
# Search a new motif
my_parser.add_argument("-m", "--motif", nargs='+',
                    type=str,
                    help="Motif to search and minimal length. ",
                    required=False)

# Execute the parse_args() method
args = my_parser.parse_args()

# Dictionary of single letter code of aa
aa_code = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
 'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}


def ds_bond(path, prot_name, distance=8):
    """
    This function obtains potential disulfide bonds between sulfur atoms in a protein sequence from a pdb file.
            Parameters:
                    path (str): relative path to .pdb file.
                    prot_name (str): the file protein name
                    distance (float): the maximal distance between sulfur atoms that form disulfide bonds
                                    (in Armstrongs)
            Returns:
                    prot_dic (dict): dictionary with the protein name, number of disulfide bonds in
                                    the whole protein, and the bonds themselves as lists with the
                                    implicated Cys residues IDs, the distance between their sulfur
                                    atoms, and the chain and model they belong to.
    """
    # Lists of bonds
    bonds=[]
    # Number of bonds
    c=0
    # Create the PDB parser and ignore warnings
    parser = PDB.PDBParser(QUIET=True)
    # Get the structure from the file
    struct=parser.get_structure('protein',path)
    # For each model, get each chain
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
                        # Distance between sulfur atoms less than 8 A for unique pairs of Cys
                        if (cys_1['SG']-cys_2['SG'])<distance and [cys_1.get_id()[1], cys_2.get_id()[1]] \
                                not in pairs and [cys_2.get_id()[1], cys_1.get_id()[1]] not in pairs:
                            pairs.append([cys_1.get_id()[1], cys_2.get_id()[1]])
                            bonds.append([cys_1.get_id()[1], cys_2.get_id()[1], cys_1['SG']-cys_2['SG'],
                                          model, chain])
                            c=c+1
    # Output dictionary
    prot_dic={'name':prot_name, 'num_bonds': c, 'di_bonds':bonds}
    print(prot_dic)



def al_helix(path, prot_name, reg_exp="\B[HKR].{0,4}[HKR]{1}[^HRKPSG]{1,3}[DE]{1}.{0,4}[DE]\B"):
    """
    This function obtains potential alpha helices in a protein from a pdb file according to a pattern
            Parameters:
                    path (str): relative path to .pdb file
                    prot_name (str): the file protein name
                    reg_exp (str): regular expression that describes the sequence pattern to search in the
                                   protein sequence. By default
            Returns:
                    prot_dic (dict): dictionary with the protein name, number of alpha helices found under
                                    the pattern given
                     in
                                    the whole protein, and the bonds themselves as lists with the
                                    implicated Cys residues IDs, the distance between their sulfur
                                    atoms, and the chain and model they belong to.
    """
    c=0
    motifs = []
    parser = PDB.PDBParser(QUIET=True)
    struct = parser.get_structure('protein', path)
    # iterate each model, chain, and residue
    # printing out the sequence for each chain
    for model in struct:
        for chain in model:
            aa_seq = []
            helix_seqs = []
            for residue in chain:
                if residue.resname in aa_code.keys():
                    aa_seq.append(aa_code[residue.resname])
                else:
                    aa_seq.append('-')
            aa_seq = ''.join(aa_seq)
            aa_seq = str(aa_seq)

            matches = re.findall(reg_exp, aa_seq)
            for helix in matches:
                if len(helix) >= 10:
                    helix_seqs.append(helix)
                    c=c+1
                    # print(helix)
            motifs.append([chain, helix_seqs])
    prot_dic = {'name': prot_name, 'num_helix': c, 'helix_seqs': motifs}
    print(prot_dic)

def b_sheets(path, name, reg_exp="\B.{1,4}[V,I,T,F,Y,W]+.{1,4}[V,I,T,F,Y,W]+.{1,4}"):
    c=0
    motifs = []
    # Create the PDB parser and ignore warnings
    parser = PDB.PDBParser(QUIET=True)
    # Get the structure from the file
    struct = parser.get_structure('protein', path)
    # iterate each model, chain, and residue
    # printing out the sequence for each chain
    for model in struct:
        for chain in model:
            aa_seq = []
            sheets_seqs = []
            for residue in chain:
                if residue.resname in aa_code.keys():
                    aa_seq.append(aa_code[residue.resname])
                else:
                    aa_seq.append('-')
            aa_seq = ''.join(aa_seq)
            aa_seq = str(aa_seq)

            matches = re.findall(reg_exp, aa_seq)
            for sheet in matches:
                if len(sheet) >= 3:
                    sheets_seqs.append(sheet)
                    c=c+1
            motifs.append([chain, sheets_seqs])
    prot_dic = {'name': prot_name, 'num_sheets': c, 'sheets_seqs': motifs}
    print(prot_dic)

def others(path, name, size=10, reg_exp="\B[V,I,T,F,Y,W]+[G,T,H,O]{1,3}"):
    c=0
    motifs = []
    # Create the PDB parser and ignore warnings
    parser = PDB.PDBParser(QUIET=True)
    # Get the structure from the file
    struct = parser.get_structure('protein', path)
    # iterate each model, chain, and residue
    # printing out the sequence for each chain
    for model in struct:
        for chain in model:
            aa_seq = []
            motif_seqs = []
            for residue in chain:
                if residue.resname in aa_code.keys():
                    aa_seq.append(aa_code[residue.resname])
                else:
                    aa_seq.append('-')
            aa_seq = ''.join(aa_seq)
            aa_seq = str(aa_seq)

            matches = re.findall(reg_exp, aa_seq)
            for sheet in matches:
                if len(sheet) >= size:
                    motif_seqs.append(sheet)
                    c=c+1
            motifs.append([chain, motif_seqs])
    prot_dic = {'name': prot_name, 'num_motif': c, 'motif_seqs': motifs}
    print(prot_dic)



if not args.alpha and not args.disulfide and not args.beta and not args.motif:
    print('At least one type of motif is required')

for path in args.input:
    try:
        if path.endswith('.pdb'):
            # Parse protein name from the path file
            prot_name = str(path).split('/')[-1]
            prot_name = str(prot_name).split('.')[0]
            if not args.disulfide==None:
                if args.disulfide == -1:
                    ds_bond(path, prot_name)
                else:
                    ds_bond(path, prot_name, args.disulfide)

            if not args.alpha == None:
                if args.alpha == 'default':
                    al_helix(path, prot_name)
                else:
                    al_helix(path, prot_name, args.alpha)

            if not args.beta == None:
                if args.beta == 'default':
                    b_sheets(path, prot_name)
                else:
                    b_sheets(path, prot_name, args.beta)

            if not args.motif == None:
                if args.motif == 'default':
                    others(path, prot_name)
                else:
                    others(path, prot_name, int(args.motif[1]), args.motif[0])

        else:
            print(path,' file must have .pdb format')
    except FileNotFoundError as ex:
        print(path,' : File not found')



