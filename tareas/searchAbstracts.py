'''
NAME
    searchAbstracts.py
VERSION
    [1.0]
AUTHOR
    Daianna Gonzalez Padilla <daianna@lcg.unam.mx>
DESCRIPTION
    This programs has two funcions, the first one search articles given an author's
    name and key words of the atricle title in a database and the second one gets a
    file with IDs of articles and returns their abstracts and IDs of articles that
    cite them
CATEGORY
     NCBI data bases analysis
USAGE
    None
ARGUMENTS
    This program doesn't take arguments
INPUT
    The input depends on the fuction used 
OUTPUT
    The output depends on the function
    
EXAMPLES
                             
    Example 1.2: gets  get_abstracts() whit the following inputs:
                        name= Otto Geiger
                        word_title= lipid
                        word_title= dna
            and the first function returns an output file called search_results_Otto_Geiger.txt with the
            following text:
                        Results:
                        7 article(s) found, whose ID(s) (for the first 20) are :
                        ['27760387', '25925947', '25040623', '20153447', '20018679', '17977153', '15341653']

            the second function returns a file called search_articles_Otto_Geiger.txt with:
                
                1. Biochim Biophys Acta Mol Cell Biol Lipids. 2017 Nov;1862(11):1287-1299. doi:
                10.1016/j.bbalip.2016.10.007. Epub 2016 Oct 17.

                Bacterial lipid diversity.

                López-Lara IM(1), Geiger O(2).

                Author information: 
                (1)Centro de Ciencias Genómicas, Universidad Nacional Autónoma de México, Av.
                Universidad s/n, Apdo, Postal 565-A, Cuernavaca, Morelos CP62210, Mexico. [...]

                The glycerophospholipids phosphatidylethanolamine, phosphatidylglycerol (PG), and
                cardiolipin (CL) are [...].

                Copyright © 2016 Elsevier B.V. All rights reserved.

                DOI: 10.1016/j.bbalip.2016.10.007 
                PMID: 27760387  [Indexed for MEDLINE]


                2. J Biol Chem. 2015 Jun 12;290(24):15102-11. doi: 10.1074/jbc.M115.639575. Epub
                2015 Apr 29.

                OlsG (Sinac_1600) Is an Ornithine Lipid N-Methyltransferase from the
                Planctomycete Singulisphaera acidiphila.

                [...]

                Artículo con ID: 27760387 citado en los artículos con IDs: ['34332123', '34093456',...
                Artículo con ID: 25925947 citado en los artículos con IDs: ['33445571', '30367006',...
                Artículo con ID: 25040623 citado en los artículos con IDs: ['33750904', '33692768',...
                [...]
                                
                
SOURCE
    https://github.com/daianna21/python_class/blob/master/tareas/searchAbstracts.py
'''
from Bio import Entrez

# Annex email
Entrez.email = "daianna@lcg.unam.mx"

def search_articles (data_base):
    """
    This function gets a data base name of NCBI and ask for an author's name and two key
    words of articles titles and returns an output file with the numer of articules found
    that match the author and title words and the IDs of the first 20 ones.
    
            Parameters:
                    data_base (str): name of the database in which to search
            Returns:
                    search_results_{Author's name}.txt (file): output file with the information
                    0 (int): if there are no matches in the search
                    the output file name (str): if there at least one match
    """
    # List of the two title words
    words_title=[]

    # Ask for the name and the words 
    name=str(input('Introduce the author\'\s name: '))
    words_title.append(str(input('Introduce the title word: ')))
    words_title.append(str(input('Introduce the title word: ')))

     # Create output file in docs/ directory
    output_file=open(f"../docs/search_results_{name}.txt", "x")

    # Term to search in the database
    term_search=f'(({name}[AUTH] OR {name}[All Fields]) AND ({words_title[0]}[Title] OR {words_title[1]}[Title]))'
    handle = Entrez.esearch(db=data_base, term=term_search) 
    result = Entrez.read(handle)

    # If there are not matches then print a message and write in the output file, return 0
    if int(result["Count"])==0:
        print('No matches.')
        output_file.write('No articles found')
        return(0)
    #Get the number and IDs list of the articles found
    else:
        output_file.write('Results:\n')
        output_file.write(result["Count"])
        output_file.write(" article(s) found, whose ID(s) (for the first 20) are : ")
        output_file.write( str(result["IdList"]))
        # Return the output file name
        return(f'search_results_{name}.txt')
        
    output_file.close()
    






def get_abstracts (data_base='pubmed'):
    """
    This function calls the function search_results to get a list of IDs of articles of a certain
    author about certain topics and then searches those IDs in pucmed to extract their abstracts
    and the IDs of articles that have cited them. 
    
            Parameters:
                    data_base (str): Defined database 'pucmed' to search there
            Returns:
                    search_articles_{Author's name}.txt (file): output file with the information
    """
    # Call search_articles function
    input_file=search_articles(data_base)
    # If there are matches in the first search of articles' IDs
    if not input_file==0:
        # Input file es the one obtained from search_articles
        input_file_name = f"../docs/{input_file}"
        # Get the author's name 
        name=str(input_file).split('results_')[1]
        name=str(name).split('.txt')[0]
        name=str(name).replace(' ','_')

        # Get the text from the input file
        input_file=open(input_file_name,"r")
        all_lines = input_file.readlines()
        input_file.close()

        # Get IDs of the input file
        ids_list=[]
        all_lines=str(all_lines)
        # Separate IDs from the rest of the text
        all_lines=(all_lines.split(': ['))
        ids=str(all_lines[1]).split(',')
        # Last ID is followed by ]" in the text
        ids[-1]=(ids[-1].split(']"'))[0]
        # Create a string like "ID1,1D2,.." for efecth
        s='"'
        for id in ids:
            # IDs are in the list as 'ID1', 'ID2',...
            # Separate from ' ' and join them with comma
            id=id.split("'")[1]
            s=str(s)+str(id)+','
            # Create a list of IDs for elink 
            ids_list.append(id)
        # String ends with "
        s=str(s)+'"'

        # Create handle to search the articles' abstracts and save them in the file out_handle 
        out_handle = open(f"../docs/search_articles_{name}.txt", "w", encoding="utf-8")
        handle = Entrez.efetch( db=data_base, id=s, rettype="abstract", retmode="text")
        data = handle.read()
        handle.close()
        out_handle.write(data)
        # Create another handle to search the ID of the articles that cited each article
        # Link name specified for pubmed
        cites=Entrez.read(Entrez.elink(dbfrom=data_base, db='pubmed', id=ids_list, LinkName='pubmed_pubmed_citedin'))
        citedin=''
        # For each article get the links of the IDs that cite it:
        # get its part of the handle, then search in LinkSetDb within the first element and search for Link
        # finally get the ID of each Link
        try:
            for i in range(0,len(ids_list)):
                citedin=f'Artículo con ID: {ids_list[i]} citado en los artículos con IDs: {[link["Id"] for link in cites[i]["LinkSetDb"][0]["Link"]]}'
                # Write the IDs in the file
                out_handle.write(citedin)
                out_handle.write('\n')
        # If there are not cites:
        except:
            print('None')
        out_handle.close()


    
    # There were not articles found for the given autor and topics    
    else:
        print('No matches.')


# Example 1:
get_abstracts()
