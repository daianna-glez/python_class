'''
NAME
    Entrez_search.py
VERSION
    [1.0]
AUTHOR
    Daianna Gonzalez Padilla <daianna@lcg.unam.mx>
DESCRIPTION
    This programs has two funcions, the first one returns the descriptions of certain elements of a given
    field for the information of a certain database. The second function search articles given an author's
    name and key words of the atricle title in a database.
CATEGORY
     NCBI data bases analysis
USAGE
    None
ARGUMENTS
    This program doesn't take arguments
INPUT
    The input depends on the fuction used. Both requiere the database name to search information
OUTPUT
    The output depends on the function
    
EXAMPLES
    Example 1.1: gets  dictionary={"FieldList":'ECNO', "LinkList":"protein_protein_small_genome"}
                       get_descriptions('protein',dictionary)
               and returns:  EC number for enzyme or CAS registry number
                             All proteins from this genome
                             
    Example 1.2: gets  search_articles(data_base="pubmed") whit the following inputs:
                        name= Valeria Mateo-Estrada
                        word_title= Genomic
                        word_title= DNA
                and returns an output file with the following text:
                        Results:
                        1 article(s) found, whose ID(s) (for the first 20) are :['34282943']
                        
                
SOURCE
    https://github.com/daianna21/python_class/blob/master/tareas/Entrez_search.py
'''



from Bio import Entrez
# Annex email
Entrez.email = "daianna@lcg.unam.mx"

def get_descriptions (data_base, fields_names):
    """
    This function gets a data base name of NCBI and a dictionary with the names of the elements
    to search and in which field of the database to search them and returns the description of
    those elements on screen. 
            Parameters:
                    data_base (str): name of the database in which to search
                    fields_names (dict): dictionary whose keys are the fields names and its
                    values the elements names
            Returns:
                    Prints the description of each element 
    """
    # Create the handle with einfo and the given database
    handle = Entrez.einfo(db = data_base)
    # Read the handle through record
    record = Entrez.read(handle)

    # For each field (key) of the dictionary get the DbInfo     
    for key in fields_names.keys():
        for field in record["DbInfo"][key]:
            # Get the description of the element whose name is the key value
            if(field["Name"]==fields_names[key]):
                print(field['Description'])
         
    
# Example 1.1:
dictionary={"FieldList":'ECNO', "LinkList":"protein_protein_small_genome"}
get_descriptions('protein',dictionary)



def search_articles (data_base):
    """
    This function gets a data base name of NCBI and ask for an author's name and two key
    words of articles titles and returns an output file with the numer of articules found
    that match the author and title words and the IDs of the first 20 ones.
    
            Parameters:
                    data_base (str): name of the database in which to search
            Returns:
                    search_results.txt (file): output file with the information
    """
    # List of the two title words
    words_title=[]
    # Create output file in docs/ directory
    output_file=open(f"../docs/search_results.txt", "x")

    # Ask for the name and the words 
    name=str(input('Introduce the author\'\s name: '))
    words_title.append(str(input('Introduce the title word: ')))
    words_title.append(str(input('Introduce the title word: ')))

    # Term to search in the database
    term_search=f'(({name}[AUTH] OR {name}[All Fields]) AND ({words_title[0]}[Title] OR {words_title[1]}[Title]))'
    handle = Entrez.esearch(db=data_base, term=term_search) 
    result = Entrez.read(handle)

    # If there are not matches then print a message and write in the output file
    if int(result["Count"])==0:
        print('No matches.')
        output_file.write('No articles found')
    #Get the number and IDs list of the articles found
    else:
        output_file.write('Results:\n')
        output_file.write(result["Count"])
        output_file.write(" article(s) found, whose ID(s) (for the first 20) are : ")
        output_file.write( str(result["IdList"]))
        
    output_file.close()

# Example 1.2:
search_articles(data_base="pubmed")

