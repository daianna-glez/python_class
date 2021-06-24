'''
NAME
    Translation.py

VERSION
    [1.0]

AUTHOR
    Daianna Gonzalez Padilla <daianna@lcg.unam.mx>

DESCRIPTION
    This programs gets a dna sequence and returns the peptide its ORF codifies if it has one

CATEGORY
     DNA sequence analysis

USAGE
    None

ARGUMENTS
    None

INPUT
    The dna sequence given by the user with only one ORF

OUTPUT
    The sequence of amino acids of the peptide the dna sequence codifies

EXAMPLES
    Example 1: gets AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA and returns MAMAPRTEINSTRING

GITHUB LINK
    https://github.com/daianna21/python_class/blob/master/scripts/Translation.py

'''
import re

# List of amino acids of the dna sequence
protein = []
# Dictionary with the codons and their corresponing amino acid: stop codons are empty since they don't carry any amino acid
gencode = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M', 'ACA': 'T',
    'ACC': 'T', 'ACG': 'T', 'ACT': 'T', 'AAC': 'N', 'AAT': 'N',
    'AAA': 'K', 'AAG': 'K', 'AGC': 'S', 'AGT': 'S', 'AGA': 'R',
    'AGG': 'R', 'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P', 'CAC': 'H',
    'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGA': 'R', 'CGC': 'R',
    'CGG': 'R', 'CGT': 'R', 'GTA': 'V', 'GTC': 'V', 'GTG': 'V',
    'GTT': 'V', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E', 'GGA': 'G',
    'GGC': 'G', 'GGG': 'G', 'GGT': 'G', 'TCA': 'S', 'TCC': 'S',
    'TCG': 'S', 'TCT': 'S', 'TTC': 'F', 'TTT': 'F', 'TTA': 'L',
    'TTG': 'L', 'TAC': 'Y', 'TAT': 'Y', 'TAA': '', 'TAG': '',
    'TGC': 'C', 'TGT': 'C', 'TGA': '', 'TGG': 'W'}

dna = input("Insert the dna sequence, make sure there is only one ORF within it:\n")
dna = dna.upper()
# Find the orf within the sequence
orf = re.search(r"ATG[ATGC]+(TAA|TGA|TAG)", dna)
# If there's a match (possible ORF) and only A,T,C,G nucleotides, the codons are obtained
if orf and not re.search(r"[^ATGC]", dna):
    # Now dna is the orf and not the whole sequence
    dna = orf.group()
    # Only if the orf is divided by 3, then their codons are translated
    if len(dna) % 3 == 0:
        # Get the codons from the start codon to the stop codon
        codons = [dna[i:i + 3] for i in range(0, len(dna), 3)]
        # The amino acid of each codon is obtained and stored in protein
        for codon_dna in codons:
            codon_dna = gencode.get(codon_dna)
            protein.append(codon_dna)
        # The elements of protein are joined to create a single string
        protein = ''.join(protein)
        print("The pepetide of the ORF within your dna sequence is: ", protein)
    else:
        print("There's no ORF within your sequence.")


else:
<<<<<<< HEAD
    print("There's no ORF within your sequence or your sequence is not dna.")
=======
    print("There's no ORF within your sequence or your sequence is not dna.")
>>>>>>> 673391c09103b048ceb5f0a21f6b1b7a2619a303
