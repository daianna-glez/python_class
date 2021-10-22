'''
NAME
    searchGOs.py
VERSION
    [1.0]
AUTHOR
    Daianna Gonzalez Padilla <daianna@lcg.unam.mx>
DESCRIPTION
    This programs receives a list of UniProt ID's and a list of GO terms and returns an output
    file with the SwissProt infomation of those ID's, including theis associated GO's if they
    exist, and their PROSITE IDs and documentation.
CATEGORY
     SwissProt and PROSITE data bases analysis
USAGE
    None
ARGUMENTS
    This program doesn't take arguments
INPUT
    List of UniProt ID's and a list of GO terms.
OUTPUT
    File with data about the ID's proteins
    
EXAMPLES
                             
    Example 1.1: gets  IDs= ["A0A0K2RVI7_9BETC", "A8R4D4_9BETC",
                             "POLG_YEFV1", "POLG_DEN1W",
                             "Q6W352_CVEN9", "D9SV67_CLOC7",
                             "A9KSF7_LACP7", "B8I7R6_RUMCH"]

                    GO_terms = ["GO:0046755", "GO:0046761",
                                "GO:0046760", "GO:0039702",
                                "GO:0046765", "GO:0046762"]

                    search_GOs(IDs,GO_terms)
                    
            and returns a file called search_GOs.txt with the following data:
            
                    ID: A0A0K2RVI7_9BETC
                    Protein name: Envelope small membrane protein 
                    GO:0046760, meaning: P:viral budding from Golgi membrane
                    Organism: Equine coronavirus.
                    SUBCELLULAR LOCATION: Host Golgi apparatus membrane 
                    Abstracts of references: 
                    No references found.
                    Prosite documentation: 
                    No prosite IDs found.

                    [...]

                    ID: POLG_YEFV1
                    Protein name: Genome polyprotein
                    GO:0046762, meaning: P:viral budding from endoplasmic reticulum membrane
                    Organism: Yellow fever virus (strain 17D vaccine) (YFV).
                    SUBCELLULAR LOCATION: [Capsid protein C]: Virion 
                    SUBCELLULAR LOCATION: [Peptide pr]: Secreted 
                    [...] 
                    Abstracts of references: 

                    1. Virus Res. 1995 Jan;35(1):35-41.

                    Complete nucleotide sequence of yellow fever virus vaccine strains 17DD and
                    17D-213.
                    dos Santos CN(1), Post PR, Carvalho R, Ferreira II, Rice CM, Galler R.
                    [...]

                    Prosite documentation: 
                    Prodoc ID: PDOC51527
                    Documentation: 
                    *************************************************************
                    * Flavivirus NS2B and NS3 protease (NS3pro) domain profiles *
                    *************************************************************
                    [...]


SOURCE
    https://github.com/daianna21/python_class/blob/master/tareas/searchGOs.py
'''

from Bio import Entrez, SeqIO
from Bio import ExPASy
from Bio.ExPASy import Prosite
from Bio.ExPASy import Prodoc
from Bio import SwissProt
import re

# User for NCBI
Entrez.email = "daianna@lcg.unam.mx" 


def search_GOs(IDs,GO_terms):
    """
    This function gets list of ID's and GO terms to search them in SwissProt and
    returns a series of interest information about those proteins. 
    
            Parameters:
                    IDs(list of str): List of UniProt ID's 
                    GO_terms(list of str): list of GO terms in the format 'GO:[number]'
            Returns:
                    search_articles_GOs.txt (file): output file with each ID, protein name,
                    GOs and their meaning, organism, subcellular locations, 3 abstracts of
                    articles from which information was extracted and domains/families
                    information from Prosite.
    """
    
    # Create an output file 
    file=open("../docs/search_GOs.txt", "w")
    # Create a handle for each ID with ExPASy and read them as objects with SwissProt
    for ID in IDs:
        handle = ExPASy.get_sprot_raw(ID)
        record = SwissProt.read(handle)
        # Write the ID
        file.write('ID: '+ID+'\n')
        # Search the full protein name within each description:
        # Letters and blank spaces only
        match=re.search(r'Name: Full=(\w|\s)+',record.description)
        match=record.description[match.start():match.end()]
        # Parse the string to obtain only the protein name and write it
        prot_name=match.split('=')[1]
        file.write('Protein name: '+prot_name+'\n')

        # Get the tuples in cross_references if they have a GO term and write their elements
        tuples=[tuple for tuple in record.cross_references for GO_term in GO_terms if GO_term in tuple]
        if len(tuples)>0:
            for tuple in tuples:
                file.write(tuple[1]+', meaning: '+tuple[2]+'\n')
        # If there are not GO's:
        else:
            file.write('No GOs associated.\n')
            
        # Get and write the organism
        file.write('Organism: '+record.organism+'\n')

        # Search all the subcelullar locations 
        # For each comment in comemets get the ones with the word 'SUBCELULLAR'
        sub_locs=[record.comments[i] for i in range(0,len(record.comments)) if 'SUBCELLULAR' in record.comments[i]]
        if len(sub_locs)>0:
             # For all the comments with the word, parse to get only the protein names and their locations
            for comment in sub_locs:
                sub_loc=comment.split('{')[0]
                file.write(sub_loc+'\n')
        else:
            file.write('No subcellular loctions found.\n')

        # Get the first three references of articles that mention the protein
        references=[]
        for reference in record.references[1:4]:
            # Append the IDs of the articles 
            references.append(reference.references[0][1])
        # Create a string with the IDs in the format: "ID1,ID2,..." for efetch 
        s='"'
        for id in references:
            # IDs are in the list as 'ID1', 'ID2',...
            # Separate from ' ' and join them with comma
            s=str(s)+str(id)+','
        # String ends with "
        s=str(s)+'"'

        # Create a handle to get all the articles abstracts and write them
        file.write('Abstracts of references: \n')
        try:
            handle = Entrez.efetch( db='pubmed', id=s, rettype="abstract", retmode="text")
            data = handle.read()
            handle.close()
            file.write(data+'\n')
        # If there are not references, print a message
        except:
            file.write("No references found."+'\n')

       # Get Prosite IDs for each record:
        IDs_ps=[tuple[1] for tuple in record.cross_references if 'PROSITE' in tuple]
       # Prosite documentation
        file.write('Prosite documentation: \n')
        # Write a message if there are not prosite IDs
        if len(IDs_ps)==0:
         file.write('No prosite IDs found.\n\n\n')
        else:
         record_pdocs=[]
        # For each prosite ID, get the information from Prosite and parse it
         for ID_ps in IDs_ps:
            handle = ExPASy.get_prosite_raw(ID_ps)
            record=Prosite.read(handle)
            # Get the unique record ids to access prosite documentation
            if record.pdoc not in record_pdocs:
                record_pdocs.append(record.pdoc)
         for record_pdoc in record_pdocs:
            file.write('Prodoc ID: '+record_pdoc+'\n')
            # Get the info and parse it
            handle=ExPASy.get_prosite_raw(record_pdoc)
            record = Prodoc.read(handle)
            # Add the documentation for that prosite ID in the file
            file.write('Documentation: \n'+ record.text +'\n\n')
            
        file.write('\n\n\n') 
 
    file.close()

# Example 1.1:
IDs= ["A0A0K2RVI7_9BETC", "A8R4D4_9BETC",
      "POLG_YEFV1", "POLG_DEN1W",
      "Q6W352_CVEN9", "D9SV67_CLOC7",
      "A9KSF7_LACP7", "B8I7R6_RUMCH"]

GO_terms = ["GO:0046755", "GO:0046761",
            "GO:0046760", "GO:0039702",
            "GO:0046765", "GO:0046762"]

search_GOs(IDs,GO_terms)


