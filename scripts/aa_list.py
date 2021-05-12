'''
NAME
    aa_list.py

VERSION
    [1.0]

AUTHOR
    Daianna Gonzalez Padilla <daianna@lcg.unam.mx>

DESCRIPTION
    This programs generates a function to determine the percentage of amino acids given in a list of a
    protein sequence

CATEGORY
     Protein sequence analysis

USAGE
    None

ARGUMENTS
    This program doesn't take arguments

INPUT
    The protein sequence and the list of amino acids

OUTPUT
    The percentage of those amino acids content in the given protein sequence

EXAMPLES
    Example 1: gets  = get_aa_percentage("MSRSLLLRFLLFLLLLPPLP", ["M"]) and returns
    "The content of those amino acids in your protein is 5.0%"

GITHUB LINK


'''

#Function to calculate the content of amino acids in a protein
#Receives the protein sequence, the amino acids list and the figures to round
def get_aa_percentage(protein_seq, aa_list, sig_figs=2):
    length = len(protein_seq)
    aa_content=0
    #Add the content (percentage) of each amino acid
    for aa in aa_list:
      aa_content=aa_content + ( round(protein_seq.count(aa) / length *100  , sig_figs))
    return(aa_content)

protein_seq=input("Enter the protein sequence:\n")
aa_list=input("Enter the amino acids list:\n")

#Evaluate these probes are true
assert get_aa_percentage("MSRSLLLRFLLFLLLLPPLP", ["M"]) == 5
assert get_aa_percentage("MSRSLLLRFLLFLLLLPPLP", ['M', 'L']) == 55
assert get_aa_percentage("MSRSLLLRFLLFLLLLPPLP", ['F', 'S', 'L']) == 70
assert get_aa_percentage("MSRSLLLRFLLFLLLLPPLP",['A','I','L','M','F','W','Y','V']) == 65

#Call the function with the correct arguments
print("The content of those amino acids in your protein is ",get_aa_percentage(protein_seq, aa_list),"%")