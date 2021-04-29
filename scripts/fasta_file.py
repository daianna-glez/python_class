'''
NAME
    fasta_file.py

VERSION
    [1.0]

AUTHOR
    Daianna González Padilla <daianna@lcg.unam.mx>

DESCRIPTION
    This programs gets a text file with dna sequences as input and returns a  file with the sequences
    in fasta format, without non alphanumeric characters, with capital letters and including the headers
    for each sequence

CATEGORY
     DNA sequence analysis

USAGE
    None

ARGUMENTS
    This program doesn't take arguments

INPUT
    The file with the DNA sequence given by the user

OUTPUT
    A fasta file with the sequences of dna

EXAMPLES
    Example 1: gets seq_1 = "ATCGTACGATCGATCGATCGCTAGACGTATCG" and returns
    > seq_1
    ATCGTACGATCGATCGATCGCTAGACGTATCG


'''

#Add re library to use regular expresions
import re

#Get the input file and create the output file
input_file_name = input("Insert the name of your text file:\n")
input_file=open(input_file_name,"r")
output_file=open("../docs/dna_sequences.fasta", "x")
# All lines of the input file are stored
all_lines = input_file.readlines()
# Close input file
input_file.close()

# For each line a header and the sequence are added in the output file
for line in all_lines:
  #The elements of each line are splitted, the first one is the name of each dna sequence
  line=line.split(' ')
  #The header for each line starts with '>'
  header='>'+ line[0]+'\n'
  output_file.write(header)
  #The third element of each line is the dna sequence
  #Remove all non alphanumeric characters and transform lowercase to uppercase
  line[2]=re.sub( '[\W*]','', line[2])
  line[2]=line[2].upper()
  output_file.write(line[2])
  output_file.write('\n')

#Close the output file
output_file.close()