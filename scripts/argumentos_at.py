'''
NAME
    argumentos_at.py

VERSION
    [1.0]

AUTHOR
    Daianna Gonzalez Padilla <daianna@lcg.unam.mx>

DESCRIPTION
    This programs gets a file with one or more dna sequences and returns an output file with the AT content
    of each sequence, from the command line

CATEGORY
     DNA sequence analysis

USAGE
    argumentos_at.py  -i input_file_path  -o output_file_path  -r sig_figs

ARGUMENTS
    -i, --input INPUT
                File with gene sequences
    -o, --output OUTPUT
                Path for the output file
    -r, --round ROUND
                Number of digits to round

INPUT
    The file with the DNA sequences given by the user

OUTPUT
    A file with the AT content of each sequence given in the input file

EXAMPLES
    Example 1: gets a file with seq_1 = "ATCGTACGATCGATCGATCGCTAGACGTATCG"
                                seq_2 = "actgatcgacgatcgatcgatcacgact"
                                seq_3 = "ACTGAC-ACTGT-ACTGTA----CATGTG"
                                seq_4 = "ATTCTGNNNNNNNNNNNNNGTC"
    and returns a new file with
                                AT content for seq_1  is 50.0
                                AT content for seq_2  is 50.0
                                AT content for seq_3  is 56.5217

GITHUB LINK
    https://github.com/daianna21/python_class/blob/master/scripts/argumentos_at.py

'''
import argparse
import os
import re

# Create the parser
my_parser = argparse.ArgumentParser(description="Script that calculates AT content using command line arguments")
# Add the arguments, all are necessary
# Add an argument to request the input file
my_parser.add_argument("-i", "--input",
                    type=str,
                    help="File with gene sequences",
                    required=True)
# Add an argument to save the output in a new file
my_parser.add_argument("-o", "--output",
                    help="Path for the output file",
                    required=True)
# Add an argument for sig figs and change it to numeric
my_parser.add_argument("-r", "--round",
                    help="Number of digits to round",
                    type=int,
                    required=True)


# Function to calculate AT content of a dna sequence
def get_at_content(dna, r):
    if dna.count("N") > 0:
        raise ValueError(f'Sequence contains {dna.count("N")} N\'s')
    length = len(dna)
    a_count = dna.count('A')
    t_count = dna.count('T')
    at_content = round((a_count + t_count) / length *100, r)
    return at_content

# Execute the parse_args() method
args = my_parser.parse_args()
#Define the input and output files, and the sig figs of output
input_file = args.input
output_file = args.output
r= args.round

#Function to validate the given paths for input and output files
def valid_path (input_file, output_file, r ):
    #If paths are valid, the AT content for each sequence is calculated
    try:
        sequences = open(input_file, 'r')
        AT_content = open(output_file, 'x')
        for sequence in sequences:
            #Get the dna sequence without non alphanumeric characters, all in uppercase
            seq=sequence.split("=")
            dna=re.sub( '[\W*]','', seq[1]).upper()
            #If there are not 'N's in the sequence, the AT content is printed and stored in the output file
            try:
                AT_content.write('AT content for ' + seq[0] + ' is ' + str(get_at_content(dna , r)))
                AT_content.write('\n')
                print('AT content for ' + seq[0] + ' is ' + str(get_at_content(dna, r)))
            except ValueError as ex:
                print('skipping invalid ',seq[0],':', ex.args[0])

        AT_content.close()
        sequences.close()

    #If the path does not exist, the path is requested again
    except FileNotFoundError as ex:
        print("Sorry, the path is not valid: ")
        input2 = input("Insert the input file path again:")
        valid_path(input2, output_file, r)
    #If the path given is a directory its content is printed, the path is requested again
    except IsADirectoryError as ex:
        print("Sorry, the path is a directory: " )
        print("Files in the directory:\n")
        print('\n'.join(os.listdir(input_file)))
        input_file2 = input("Insert the input file path again:")
        valid_path(input_file2, output_file, r)

    #If the output file name already exist another one is requested
    except FileExistsError as ex:
        print("Sorry, the output file name already exists:" )
        output_file2 = input("Insert the output file name path again:")
        valid_path(input_file, output_file2, r)

#Call the function with the arguments given
valid_path(args.input, args.output, args.round)



