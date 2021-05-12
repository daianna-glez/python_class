'''
NAME
    ATCG_content.py

VERSION
    [1.0]

AUTHOR
    Daianna González Padilla <daianna@lcg.unam.mx>

DESCRIPTION
    This programs gets a dna sequence as input and returns as output the content of each nucleotide
    of the sequence

CATEGORY
     DNA sequence analysis

USAGE
    None

ARGUMENTS
    This program doesn't take arguments

INPUT
    The DNA sequence given by the user

OUTPUT
    The results are printed on screen

EXAMPLES
    Example 1: gets ATTCGCCC and returns "The A content of the sequence is: 1"
                                         "The T content of the sequence is: 2"
                                         "The C content of the sequence is: 4"
                                         "The G content of the sequence is: 1"
GITHUB LINK
    https://github.com/daianna21/python_class/blob/master/scripts/ATCG_content.py
'''

# 1. The DNA sequence is requested from the user and stored on the dna variable
print("Insert the DNA sequence:\n")
dna = input()

# 2. For each nucleotide in the sequence, it is evaluated if the are valid, this is, if they are either A, T, C or G
# Tuple nucleotides contains the valid letters (nucleotides)
nucleotides=('A', 'T', 'C', 'G')
# For each character in the sequence given, we suppose the flag dna_nucleotides is 1, that indicates that the letter
# is valid, otherwise, the flag is 0 and that happens if the letter is not in the tuple nucleotides
for char in dna:
  dna_nucleotides=1
  if char not in nucleotides:
    print("DNA sequence can only contain A, T, C and G")
    dna_nucleotides=0
    break
# 3. If all the letters are valid, dna_nucleotides is 1 and the content of each nucleotide of the dna sequence
# is printed, it is obtained with dna.count()
if (dna_nucleotides==1):
    print('The A content of the sequence is: ', dna.count('A'))
    print('The C content of the sequence is: ', dna.count('C'))
    print('The G content of the sequence is: ', dna.count('G'))
    print('The T content of the sequence is: ', dna.count('T'))
