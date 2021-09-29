#Clase 4: Bases de datos Entrez ( E-utilities de NCBI)

#E utilities son la manera de interactuar con Entrez : URLs que se conectan a Entrez
# y da un archivo fasta, xml

from Bio import Entrez
from pprint import pprint  # para mejor visualización de diccionarios!!
# Anexar un correo
Entrez.email = "daianna@lcg.unam.mx"  # IMPORTANTE!!!

#Einfo
# Objeto handle tiene a la einfo
handle = Entrez.einfo()
#result = handle.read()
#Formato amigable para visualización
record = Entrez.read(handle)
print(record)
#chequemos qué hay en einfo
#print(result)

# obtenemos diccionario (llave "Dblist")
print(record["DbList"][0:3])  # primeras 3 bases de datos


handle.close()


#Podemos checar más información sobre las bases de datos en einfo:
handle = Entrez.einfo(db = "pubmed") # indicar db de interes
record = Entrez.read(handle)
handle.close() #cerramos handle
print(record["DbInfo"]["Description"])  # descripcion de pubmed

# para saber qué podemos consultar en DbInfo
print(record["DbInfo"].keys())
print(record["DbInfo"]["LastUpdate"])
print(record["DbInfo"]["FieldList"])



# imprimir todos los campos a los que podemos accesar de pubmed 
for field in record["DbInfo"]["FieldList"]:
  #Hay muchos fields para el Field list
  print("%(Name)s, %(FullName)s, %(Description)s" % field)
  #print(field['Name'], field['FullName'], field['Description'])




# URL del request:
## el URL contiene tool=biopython y email=cgil@lcg.unam.mx, que son requeridos y que son añadidos  por biopython.
print(handle.url)


# Ejercicio 1:
handle = Entrez.einfo(db = "genome")
record = Entrez.read(handle)
print(record["DbInfo"]["Description"])
for field in record['DbInfo']["FieldList"][0:2]:
    print(" %(Description)s" % field)
    #print(field['Description'])






# Esearch: Entrez.esearch( base de datos a buscar , termino )
## Buscaremos el termino "biopython" en PubMed, checaremos cuantos resultados obtenemos con record["Count"]
handle = Entrez.esearch(db = "pubmed", term = "biopython")
record = Entrez.read(handle) 
print(record["Count"])
handle.close()

#retmax
#Parámetro que indica número máximo de retrieves (no. de resultados que muestra del total).Modificar retmax:

print(len(record["IdList"]))  #chequemos tamaño 
count = int(record["Count"]) #cambiemos retmax por long de Counts
handle = Entrez.esearch(db="pubmed", term="biopython", retmax=count)
record = Entrez.read(handle) 
handle.close()
print(len(record["IdList"])) # ahora es de 35!!

print(record.keys())  # info que podemos obtener

#Ejemplo: buscar un nombre en el campo Autor
handle = Entrez.esearch(db="pubmed", term='Valeria Mateo-Estrada',field='AUTH')
record = Entrez.read(handle)
print(handle.url)  # URL de request
handle.close()
print(record["IdList"])  # ids de artículos
## Si queremos buscar en varios campos, podemos incluirlo todo en term usando corchetes para el campo
handle = Entrez.esearch(db="pubmed", term='Valeria Mateo-Estrada[AUTH]')




# Usando operadores booleanos: búsquedas más específicas
## Buscar combinacion de campos:

termino = "(Aedes[Title] OR Aedes[All Fields])AND((RNA-Seq[Title] OR transcriptomic[Title]) OR (transcriptome[Title] OR sequencing[Title]))"
handle = Entrez.esearch(db="pubmed", term=termino)
result = Entrez.read(handle)
print(result["Count"])  #cuantos encontró





#Tarea:
##Primera parte
##Empleando Entrez.einfo y Entrez.read, imprime la descripción de los siguientes campos de la base de datos "protein":

from Bio import Entrez
Entrez.email = "daianna@lcg.unam.mx"
handle = Entrez.einfo(db = "protein")
record = Entrez.read(handle)
for field in record["DbInfo"]["FieldList"]:
    if(field["Name"]=='ECNO'):
        print(field['Description'])
 
for link in record["DbInfo"]["LinkList"]:
    if(link["Name"]=="protein_protein_small_genome"):
        print(link['Description'])


## Segunda parte
def search_db (data_base):
    words_title=[]
    output_file=open(f"../docs/search_results.txt", "x")
    
    name=str(input('Introduce the author\'\s name: '))
    words_title.append(str(input('Introduce the title word: ')))
    words_title.append(str(input('Introduce the title word: ')))
    
    term_search=f'(({name}[AUTH] OR {name}[All Fields]) AND ({words_title[0]}[Title] OR {words_title[1]}[Title]))'
    handle = Entrez.esearch(db=data_base, term=term_search) 
    result = Entrez.read(handle)
    if int(result["Count"])==0:
        print('No matches.')
        output_file.write('No articles found')
    else:
        output_file.write(('Results:\n',result["Count"]," articles found, whose IDs are: ", result["IdList"]))
    output_file.close()

search_db(data_base="pubmed")




