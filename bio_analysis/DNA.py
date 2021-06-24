'''
NAME
    DNA.py

VERSION
    [1.0]

AUTHOR
    Daianna Gonzalez Padilla <daianna@lcg.unam.mx>

DESCRIPTION
    This module is part of the bio_analysis package.
    Gets a dna sequence and works with it through a series of functions

CATEGORY
     DNA sequence analysis

FUNCTIONS
    get_sequence(str): validates the sequence given is DNA
    get orf(str): obtains the orf of the sequence
    get_protein(str): obtains the peptide codified by the sequence
    atcg_content(str): determines the A,T,C ang G content of the sequence
    at_content(str): determines the AT content of the sequence
    gc_content(str): determines the GC content of the sequence
    at_regions (str): obtains the AT rich regions of the sequence
    revcomp_seq(str): obtains the complementary reverse sequence

GITHUB LINK
    https://github.com/daianna21/python_class/blob/master/pgr/DNA.py...

'''
import re

def get_sequence (dna):
    '''
    This function evaluates if the sequence given has only A,T,C, and G characters, the dna nucleotides
            Parameters:
                    dna (str): the sequence given by the user
            Returns:
                    dna (str): if it has only A,T,C and G
                    0 (int): otherwise
    '''
    # The sequence must be uppercase
    dna = dna.upper()
    # Look up the characters that are not A,T,G or C in the sequence
    matches = re.finditer(r"[^ATGC]", dna)
    #If all characters are A,C,T, or G, the dna sequence is returned,
    # otherwise a cero is returned and the ambiguous bases are printed
    if not re.search(r"[^ATGC]", dna):
        return(dna)
    else:
        print("Your sequence is not DNA:\n")
        for m in matches:
            base = m.group()
            pos = m.start()
            print(base, " found at position ", pos)
        exit()

def get_orf(dna):
    """
    This function obtains the ORF within a dna sequence if there is one
            Parameters:
                    dna (str): the sequence given by the user
            Returns:
                    orf.group() (str): the part of the original sequence corresponding to an ORF
                    0 (int): if there is no ORF within the sequence
    """
    #Call the function to determine if the sequence is dna
    dna = get_sequence(dna)
    # Find the orf within the sequence starting with a start codon and ending with a stop codon
    orf = re.search(r"ATG[ATGC]+(TAA|TGA|TAG)", dna)
    # If there's a match (possible ORF) and the total length of the match is divided by 3
    #  the codons are obtained
    if orf and len(orf.group()) % 3 == 0:
        return(orf.group())
    else:
        print("There is no ORF within your sequence.")
        return(0)


def get_protein(dna):
    """
    This function translates an ORF of a dna sequence into the corresponding  peptide
            Parameters:
                    dna (str): the sequence given by the user
            Returns:
                    protein (str): the resulting string of amino acids
                    0 (int): if there is no ORF within the sequence
    """
    # List of amino acids of the dna sequence
    protein = []
    # Dictionary with the codons and their corresponing amino acid: stop codons are empty since they
    # don't carry any amino acid
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
    # Call the function to determine if the sequence is dna
    dna=get_sequence(dna)
    # Find the orf within the sequence by calling the function get_orf(), if it does not return a 0, the codons are
    # obtained
    if not get_orf(dna)==0:
        # The orf is obtained
        orf = get_orf(dna)
        # Get the codons from the start codon to the stop codon
        codons = [orf[i:i + 3] for i in range(0, len(orf), 3)]
        # The amino acid of each codon is obtained and stored in protein
        for codon_dna in codons:
            codon_dna = gencode.get(codon_dna)
            protein.append(codon_dna)
        # The elements of protein are joined to create a single string
        protein = ''.join(protein)
        return (protein)
    else:
        return (0)


def atcg_content(dna):
    """
    This function calculates the content of each nucleotide of a dna sequence
            Parameters:
                    dna (str): the sequence given by the user
            Returns:
                    contents (dict): a dictionary with nucleotides as keys and their percentages as values
    """
    # Call the function to determine if the sequence is dna
    dna = get_sequence(dna)
    contents={}
    # For each nucleotide the content is determined in percentage
    for char in ['A','T','C','G']:
        contents[char]=dna.count(char)/len(dna)*100
    return (contents)


def at_content(dna):
    """
    This function calculates the AT content of a dna sequence
            Parameters:
                    dna (str): the sequence given by the user
            Returns:
                    at_content (float): the percentage of AT the sequence has
    """
    # Call the function to determine if the sequence is dna
    dna = get_sequence(dna)
    # Calculate how many A's and T's are in the sequence
    a_count = dna.count('A')
    t_count = dna.count('T')
    # Final value
    at_content = (a_count + t_count) / len(dna)*100
    return (at_content)


def gc_content(dna):
    """
    This function calculates the GC content of a dna sequence
            Parameters:
                    dna (str): the sequence given by the user
            Returns:
                    gc_content (float): the percentage of GC the sequence has
    """
    # Call the function to determine if the sequence is dna
    dna = get_sequence(dna)
    # Calculate how many G's and C's are in the sequence
    g_count = dna.count('G')
    c_count = dna.count('C')
    # Final value
    gc_content = (g_count + c_count) / len(dna)*100
    return (gc_content)

def at_regions (dna):
    """
    This function obtains the AT rich regions of a dna sequence, if there are
            Parameters:
                    dna (str): the sequence given by the user
            Returns:
                    at_regions (list): a list with all AT rich regions of the sequence
                    0 (int): if there are no AT rich regions in the sequence
    """
    at_regions=[]
    # Call the function to determine if the sequence is dna
    dna = get_sequence(dna)
    # Find the regions with 5 or more A's or T's
    regions = re.finditer(r'((A|T){5,})', dna)
    # The regions found are stored in the list
    for region in regions:
        at_regions.append(region.group())
    # If there are not AT rich regions, the list is empty and a message is printed
    if at_regions==[]:
        print("No AT rich regions found")
        return(0)
    else:
        print("AT rich regions:\n")
        return (at_regions)

def revcomp_seq(dna):
    """
    This function obtains the complementary reverse sequence of a given dna sequence
            Parameters:
                    dna (str): the sequence given by the user
            Returns:
                    com_dna.upper() (str): the resulting sequence
    """
    # Call the function to determine if the sequence is dna
    dna = get_sequence(dna)
    # The reverse sequence of the original one
    rev_dna=dna[::-1]
    # Each nucleotide is replaced by its complementary
    com_dna=rev_dna.replace('A','t').replace('T','a').replace('G','c').replace('C','g')
    return(com_dna.upper())
