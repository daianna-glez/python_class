'''
NAME
    GenBank.py
VERSION
    [1.0]
AUTHOR
    Daianna Gonzalez Padilla <daianna@lcg.unam.mx>
DESCRIPTION
    This programs takes a DNA sequence as a class and works with it through its methods and subclasses
CATEGORY
     GenBank files analysis
USAGE
    None
ARGUMENTS
    This program doesn't take arguments
INPUT
    The input of the function is the path to the file and a list of genes to find 
OUTPUT
    Returns a series of interest traits of the organisms
    
EXAMPLES
    Example 1: gets   print(info_organism("C:/Users/hp/Downloads/virus.gb", [['L'],['G']]))
               and returns  Organism:  Isfahan virus
                            Date:  13-AUG-2018
                            County:  ['Iran:Isfahan province']
                            Isolation source:  ['Phlebotomus papatasi']
                            Gene name:  ['G']
                            Protein name:  ['glycoprotein']
                            DNA:  ATGACTTCAGTCTTA ...
                            RNA:  AUGACUUCAGUCUUA ...
                            Protein:  MTSVL ...
                            Gene name:  ['L']
                            Protein name:  ['large polymerase protein']
                            DNA:  ATGGATGAGTACTCT ...
                            RNA:  AUGGAUGAGUACUCU ...
                            Protein:  MDEYS ...
SOURCE
    https://github.com/daianna21/python_class/blob/master/tareas/GenBank.py
'''

from Bio import SeqIO

def info_organism(path, genes):
    """
    This function gets a genbank file and returns a series of interest traits
    of the organisms
            Parameters:
                    path: absolute path to file
                    genes: list with the names of the genes, each one in capital
                           letters and within brackets 
            Returns:
                    Organism: name of the organism
                    Date: date whern it was reported
                    Country: country where it was reported
                    Isolation source: source of where it was isolated 
                    Gene name: single letter name of the gene
                    Protein names: genes products
                    DNA: the first 15 nucleotides of each gene
                    RNA: the transcrit of DNA
                    Protein: peptide from DNA
    """
    #Get the file and format
    for gb_record in SeqIO.parse(path, "genbank"):
        #Get metadata
        print('Organism: ',gb_record.annotations['organism'])
        print('Date: ',gb_record.annotations['date'])
        #Get country from source 
        print('County: ',gb_record.features[0].qualifiers['country'])
        print('Isolation source: ',gb_record.features[0].qualifiers['isolation_source'])

        #Get those features that are CDS
        for feature in gb_record.features:
            if feature.type=='CDS':
                #Check if the CDS correspond to one of the given genes
                for gene in genes:
                     if feature.qualifiers['gene']==gene:
                         #Get info from the feature
                         print ('Gene name: ',feature.qualifiers['gene'])
                         print ('Protein name: ',feature.qualifiers['product'])
                         start = feature.location.nofuzzy_start
                         end = start + 15
                         seq = gb_record.seq[start:end]
                         print('DNA: ',seq, '...')
                         print('RNA: ',seq.transcribe(),'...')
                         print('Protein: ',seq.translate(),'...')
                     
       
#Example
print(info_organism("C:/Users/hp/Downloads/virus.gb", [['L'],['G']]))
