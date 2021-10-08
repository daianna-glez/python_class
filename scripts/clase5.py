#Clase 5:

from Bio import Entrez

#Ejercicio 2: Linajes
#PRIMERA PARTE: esearch para buscar 1er organismo en taxonomy
handle = Entrez.esearch(db="Taxonomy", term="Notoryctes typhlops")
record = Entrez.read(handle)
record["IdList"] # obtenemos ID de taxonomia
id_taxo = record["IdList"][0] #guarda ID

# SEGUNDA PARTE: efetch para obtener archivo de taxonomia
handle = Entrez.efetch(db="Taxonomy", id=id_taxo, retmode="xml")
Notoryctes = Entrez.read(handle)
Notoryctes[0].keys()  #checamos qué informacion tenemos

# PRIMERA PARTE
handle = Entrez.esearch(db="Taxonomy", term="Chrysochloris asiatica")
record = Entrez.read(handle)
id_taxo = record["IdList"][0]
print(Notoryctes[0]["Lineage"])  #checamos linaje
#SEGUNDA PARTE
handle = Entrez.efetch(db="Taxonomy", id=id_taxo, retmode="xml")
Chryso = Entrez.read(handle)
print(Chryso[0]["Lineage"])

## Para comparar linajes: con zip que hace tuplas con los elementos 1-1:
l1=Chryso[0]["Lineage"].split(';')
l2=Notoryctes[0]["Lineage"].split(';')

def comparar (org1, org2):
    i=-1
    for lin1, lin2 in zip(org1,org2):
        i=i+1
        if not lin1==lin2:
            print(i,lin1,lin2)
            break




# Elink: buscar en otras db, permite pasar de una db a otra
ids = "15718680" # id a buscar
# elink buscara los ids de protein en la base de datos de gene
record = Entrez.read(Entrez.elink(dbfrom="protein", id=ids,db='gene'))
print(record[0]) #visualicemos record

ids = "15718680,157427902" # ids a buscar
# elink buscara los ids de protein en la base de datos de gene
record = Entrez.read(Entrez.elink(dbfrom="protein", id=ids,db='gene'))
print(record[0]) #visualicemos record


##Modificaremos URL para que IDs correspondan a los solicitados:
# Función para generar la URL como en la documentación de ENTREZ
from urllib.request import urlopen
from urllib.parse import urlencode
def elink_multiple(dbfrom, ids, db,
                  mirror="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"):
    # diccionario con lo que tendrá el URL
    parameters = {"dbfrom": dbfrom, "db":db, "id": ids, "tool":"biopython", "email":Entrez.email}
    # Creamos la URL
    command = urlencode(parameters, doseq=True)
    url = "%s/elink.fcgi?%s" % (mirror, command)
    handle = urlopen(url)
    return(handle)
pmids = ["15718680","157427902"] # ids a buscar
handle = elink_multiple(dbfrom="protein", ids=pmids, db="gene")
handle.url  # chequemos URL
record = Entrez.read(handle)
handle.close()
#print(record[1])


#Obtener citas:
pmid = "32703847" #pubmed id
results = Entrez.read(Entrez.elink(dbfrom="pubmed", db="pmc", id=pmid))

print(results)
#Nos interesa la parte de pubmed_pmc_refs podemos espeficarlo desde que hacemos request de elink:
results = Entrez.read(Entrez.elink(dbfrom="pubmed", db="pmc",
                                  LinkName="pubmed_pmc_refs", from_uid=pmid))

# Guardemos links de PMC
pmc_ids = [link["Id"] for link in results[0]["LinkSetDb"][0]["Link"]]
pmc_ids
#ahora partimos de PMC a pubmed
handle=Entrez.elink(dbfrom="pmc", db="pubmed",
                                    LinkName="pmc_refs_pubmed",
                                    id=",".join(pmc_ids))
results2 = Entrez.read(handle)
#guardamos links
pubmed_ids = [link["Id"] for link in results2[0]["LinkSetDb"][0]["Link"]]
pubmed_ids
comparar(l1,l2)
print(handle.url)



input_file_name = input("../docs/search_results.txt")
input_file=open(input_file_name,"r")
all_lines = input_file.readlines()
input_file.close()
print(all_lines)

















