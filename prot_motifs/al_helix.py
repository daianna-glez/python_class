from Bio import PDB
import argparse
import re

aa_code = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
 'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

# Amino acids
aa={'A','V','L','I','F','W','M','P','G','Y','S','T','C','N','Q','K','R','H','D','E'}
aliphatic={'I','L','V','A','G'}
aromatic={'F','Y','W'}
hidroxi_aa={'S','T'}
tio_aa={'C','M'}
positive={'H','K','R'}
negative={'D','E'}
charged=positive|negative
polar={'T','G','S','Q','Y','W','N'}|charged
apolar=aa^polar
smallest={'A','G','C','S'}
small={'P','V','T','N','D'}|smallest
hydrophobic={'M','K','T','H','F','Y','W','V','A','L','I'}


def al_helix(paths, reg_exp="\B[HKR].{0,4}[HKR]{1}[^HRKPSG]{1,3}[DE]{1}.{0,4}[DE]\B"):
    # Get each protein/file
    for path in paths:
        # Parse protein name from the path file
        prot_name = str(path).split('/')[2]
        prot_name = str(prot_name).split('.')[0]
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
                helix_seqs = []
                for residue in chain:
                    if residue.resname in aa_code.keys():
                        aa_seq.append(aa_code[residue.resname])
                    else:
                        aa_seq.append('-')
                aa_seq = ''.join(aa_seq)
                aa_seq = str(aa_seq)
                print(aa_seq)

                matches = re.findall(reg_exp, aa_seq)
                for helix in matches:
                    if len(helix) >= 10:
                        helix_seqs.append(helix)
                        c=c+1
                        # print(helix)
                motifs.append([chain, helix_seqs])
        prot_dic = {'name': prot_name, 'num_bonds': c, 'helix_seqs': motifs}
        print(prot_dic)


#r"\B[HKR].{0-3}([HKR]{1}[^HRK]{1,3}[DE]{1}){1,3}.{0,3}[DE]\B"
al_helix(['../docs/1kcw.pdb', '../docs/1fat.pdb'],r"\B[HKR].{0,4}.*")

