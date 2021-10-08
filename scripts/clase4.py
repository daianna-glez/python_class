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


# Esearch busca los IDs y al meterlos en efetch da cosas de esos IDs
# EGquery Muestra en cuales de las bases de datos podemos encontrar información de nuestra búsqueda.

termino = "(Aedes[Title] OR Aedes[All Fields])AND((RNA-Seq[Title] OR transcriptomic[Title]) OR (transcriptome[Title] OR sequencing[Title]))"
handle = Entrez.egquery(term=termino)
record = Entrez.read(handle)
for row in record["eGQueryResult"]:
    print(row["DbName"], row["Count"])


#Espell: Ayuda a corregir búsqueda con sugerencias de ortografía:
handle = Entrez.espell(term="biopythooon")
record = Entrez.read(handle)
record["Query"] # lo que añadimos
record["CorrectedQuery"]  # la sugerencia

#Esummary:
handle = Entrez.esummary(db="taxonomy", id="9913,30521")
record = Entrez.read(handle)
len(record)
record[0].keys()
record[0]["Id"]

#Efetch: Regresa records en formato especificado (tipo y modo). En esta tabla podemos ver las bases de datos en las que Efetch puede
# interactuar y sus valores para retmode y rettype

from Bio import Entrez, SeqIO  
# db = "nuccore" tambien es valido
handle = Entrez.efetch(db="nucleotide", id="HE805982", rettype="gb", retmode="text")
# leemos archivo genebank
record = SeqIO.read(handle, "genbank")
handle.close()
print(record)  # imprimimos archivo


# Guardar info en un archivo
filename = "HE805982.gb"  #nombre del archivo a generar
with Entrez.efetch(db="nucleotide",id="HE805982",rettype="gb", retmode="text") as file:
    with open(filename, "w") as handle:
        handle.write(file.read())  #escribimos archivo
# parseamos archivo con SeqIO, indiicamos que es tipo genbank
record = SeqIO.read("HE805982.gb", "genbank") 
record



# Archivos fasta con handle.read
out_handle = open("files/prueba.fasta", "w")
fetch_handle = Entrez.efetch(db="nucleotide", id="1919569438, 1919569357, 1251949171",
                            rettype="fasta", retmode="text")
data = fetch_handle.read()  #usar handle.read()
fetch_handle.close() #cerrar handle
out_handle.write(data) #escribir archivo
out_handle.close() #cerrar archivo




