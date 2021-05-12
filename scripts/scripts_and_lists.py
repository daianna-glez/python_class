'''
NAME
    scripts_and_lists.py

VERSION
    1.0

AUTHOR
    Daianna González Padilla


DESCRIPTION
        From a dna sequence as input, it identifies the position of the start codon AUG of a transcript of rna
        and the end of it. It also provides the transcript itself

CATEGORY
        Genomic Sequence


'''
#The DNA sequence is requested from the user and stored in the dna variable
dna = input("Enter the DNA sequence from which you want to get the transcript:\n")
#The RNA secuence is obtained from the DNA by replacing T with U
rna=dna.replace('T','U')
#The variable start stores the numeric value of the position at which the AUG codon is located
start=rna.find('AUG')
#The position of this codon is printed (actually the position is start+1)
print('AUG codon is at position: ',start)
#These variables (end) store the start position of the posible stop codons, +1 is added to adjust the position since
# the index starts from cero
end1=rna.find('UAA')+1
end2=rna.find('UGA')+1
end3=rna.find('UAG')+1
#The resultes are printed: if the stop codon wasn't found it is indicated, otherwise, the position and the
# transcript are given
if (end1==0):
    print("UAA stop codon is not in the transcript ")
else:
    print('UAA codon starts at position ', end1, ' and ends at position ', end1 + 2)
    print('The transcript from the dna sequence with this stop codon is: ', dna[start:end1 + 2])

if (end2==0):
    print("UGA stop codon is not in the transcript ")
else:
    print('UGA codon starts at position ', end2, ' and ends at position ', end2 + 2)
    print('The transcript from the dna sequence with this stop codon is: ', dna[start:end2 + 2])

if (end3==0):
    print("UAG stop codon is not in the transcript ")
else:
    print('UAG codon starts at position ', end3, ' and ends at position ', end3 + 2)
    print('The transcript from the dna sequence with this stop codon is: ', dna[start:end3 + 2])




