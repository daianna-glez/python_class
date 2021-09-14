'''
NAME
    POO.py

VERSION
    [1.0]

AUTHOR
    Daianna Gonzalez Padilla <daianna@lcg.unam.mx>

DESCRIPTION
    This programs takes a DNA sequence as a class and works with it through its methods and subclasses

CATEGORY
     DNA analysis

USAGE
    None

ARGUMENTS
    This program doesn't take arguments

INPUT
    The input is the DNA object, just giving the sequence

OUTPUT
    Depends on the method

EXAMPLES
    Example 1: with the following code:  seq1=DNA('GGCGCGGCATGCCCCGGGTAAAGCG') and print(seq1.__dict__)
                                                                                     print(seq1.revcomp_seq())
                                                                                     print(seq1.gc_content())
                                                                                     print(seq1.get_orf())
    and returns {'seq': 'GGCGCGGCATGCCCGGGTAAAGCG'}
                TTACCCGGGCAT
                75.0
                ATGCCCGGGTAA

SOURCE
    https://github.com/daianna21/python_class/blob/master/tareas/POO.py


'''

import re

# Create class of DNA sequence
class DNA():
    def __init__(self, seq):
        # Instance attributes
        self.seq = seq

    # Methods
    def gc_content(self):
        """
        This function calculates the GC content of a dna sequence
                Parameters:
                        seq (str): the sequence given by the user
                Returns:
                        gc_content (float): the percentage of GC the sequence has
        """
        # Calculate how many G's and C's are in the sequence
        g_count = self.seq.count('G')
        c_count = self.seq.count('C')
        # Final value
        gc_cont = (g_count + c_count) / len(self.seq) * 100
        return gc_cont

    def revcomp_seq(self):
        """
        This function obtains the complementary reverse sequence of a given dna sequence
                Parameters:
                        seq (str): the sequence given by the user
                Returns:
                        com_dna.upper() (str): the resulting sequence
        """
        # The reverse sequence of the original one
        rev_dna = self.seq[::-1]
        # Each nucleotide is replaced by its complementary
        com_dna = rev_dna.replace('A', 't').replace('T', 'a').replace('G', 'c').replace('C', 'g')
        return (com_dna.upper())


# Inheritance: subclass ORF
# This class is part of DNA but works with the ORF within the given dna sequence
class ORF(DNA):
    def get_orf(self):
        """
        This function obtains the ORF within a dna sequence if there is one
                Parameters:
                        seq (str): the sequence given by the user
                Returns:
                        orf.group() (str): the part of the original sequence corresponding to an ORF
                        0 (int): if there is no ORF within the sequence
        """
        # Find the orf within the sequence starting with a start codon and ending with a stop codon
        orf = re.search(r"ATG[ATGC]+(TAA|TGA|TAG)", self.seq)
        # If there's a match (possible ORF) and the total length of the match is divided by 3
        #  the codons are obtained
        if orf and len(orf.group()) % 3 == 0:
            return (orf.group())
        else:
            return (0)

    # Overriding and polimorfism: applying the revcomp_seq function in ORF class gives a different result than  applied
    # in DNA class
    def revcomp_seq(self):
        """
        This function obtains the complementary reverse sequence of the ORF of a given dna sequence
                Parameters:
                        seq (str): the sequence given by the user
                Returns:
                        com_dna.upper() (str): the resulting sequence
        """
        # Obtain the ORF with the previous function
        orf = self.get_orf()
        if orf == 0:
            print("There's no orf")
        else:
            # The reverse sequence of the orf
            rev_dna = orf[::-1]
            # Each nucleotide is replaced by its complementary
            com_dna = rev_dna.replace('A', 't').replace('T', 'a').replace('G', 'c').replace('C', 'g')
            return (com_dna.upper())
#Examples:
seq1=ORF('GGCGCGGCATGCCCGGGTAAAGCG')
seq2=DNA('AGTGCTGCTC')
seq3=ORF('AGTGCTGCTC')

print(seq1.__dict__)
print(seq1.revcomp_seq())
print(seq1.gc_content())
print(seq1.get_orf())


print(seq2.revcomp_seq())
print(seq2.gc_content())
print(seq3.revcomp_seq())
print(seq3.gc_content())
