'''
NAME
    AT_regions.py

VERSION
    [1.0]

AUTHOR
    Daianna Gonzalez Padilla <daianna@lcg.unam.mx>

DESCRIPTION
    This programs gets a dna sequence and returns the AT rich regions of it.

CATEGORY
     DNA sequence analysis

USAGE
    None

ARGUMENTS
    None

INPUT
    The dna sequence given by the user

OUTPUT
    Those regions of the dna sequence that are AT rich

EXAMPLES
    Example 1: gets CTGCATTATATCGTACGAAATTATACGCGCG and returns ATTATAT
                                                                AAATTATA

GITHUB LINK
    https://github.com/daianna21/python_class/blob/master/scripts/AT_regions.py

'''

import re

class AmbiguousBaseError(Exception):
    pass
def AT_rich_seq(dna):
    '''
     This function receives a dna sequence and returns the AT rich regions (the ones with 5 to 15 A or T
     nucleotides). It also evaluates that the sequence has only A,T.C or G.
    '''
    #Look up the characters that are not A,T,G or C in the sequence
    matches = re.finditer(r"[^ATGC]", dna)

    # If the sequence has only A,T,G or C, the regions that have 5 to 15 A or T are printed, otherwise
    # an error is raised and the non allowed characters and their positions are printed
    if not re.search(r"[^ATGC]", dna):
        regions = re.finditer(r'((A|T){5,15})', dna)
        empty=False
        print("AT rich regions:")
        for region in regions:
            if region.group():
                empty=True
            print(region.group())
        # If there are not AT rich regions, a message is printed
        if empty==False:
            print("No AT rich regions found")
<<<<<<< HEAD

=======
>>>>>>> 673391c09103b048ceb5f0a21f6b1b7a2619a303
    else:
        for m in matches:
            base = m.group()
            pos = m.start()
            print(base, " found at position ", pos)
        raise AmbiguousBaseError()


try:
    #Ask for the dna sequence and call the function
    dna = input("Insert the dna sequence, use capital letters:\n")
    AT_seq = AT_rich_seq(dna)
except AmbiguousBaseError:
    print('Error: Ambiguous bases found in the dna sequence')
