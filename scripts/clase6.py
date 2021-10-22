#Clase 6:
# PDB tiene estructuras 3D de rna, complejos macromoleculares por cristalografía
# o electroscopía mientras que en Uniprot no todo estña anotado por curadores
# (sin evidencia) pero tiene más info que otras db
# Prosite para ver patrones y dominios 

#Bases de datos:  ExPASy= forma de interactuar con todas las db del instituto suizo

from Bio import Entrez, SeqIO
Entrez.email = "daianna@lcg.unam.mx" # Always tell NCBI who you are
# Buscar info en esearch
# Busco mi organismo y gen en la base de datos de proteínas
handle = Entrez.esearch(db="protein", term="Aedes aegypti[Orgn] AND APY[Gene]") 
record = Entrez.read(handle) 
# ¿Cuántos resultados encontró? y el ID
print(record["Count"])
print(record['IdList'])

# Extraer info desde la base de datos de protein
handle = Entrez.efetch(db="protein", id="193806340", rettype="gb", retmode="text")
record = SeqIO.read(handle, "genbank")
handle.close()
# P50635.2 -> P50635 version 2
print(record.id)
print(record.annotations['db_source'])


# Información en otras bases de datos (ej. UniProt)
db_source = record.annotations['db_source']
apy_prot = record.annotations['accessions']

print(apy_prot)

#Con la db de gene:
# Busco la información
handle = Entrez.esearch(db="gene", term="Aedes aegypti[Orgn] AND apy[Gene]")
record = Entrez.read(handle)   
# La obtengo 
handle = Entrez.efetch(db="gene", id=record['IdList'][0], retmode="xml")
record = Entrez.read(handle, "genbank")
# Sinonimos
print(record[0]['Entrezgene_gene']['Gene-ref']['Gene-ref_syn'])
sinonimos = record[0]['Entrezgene_gene']['Gene-ref']['Gene-ref_syn']
#Sinónimo en uniprot tiene _
print([uniprot for uniprot in sinonimos if '_' in uniprot])



#Ejercicio 1:
# Busco mi organismo y gen en la base de datos de proteínas
handle = Entrez.esearch(db="protein", term="Aedes aegypti[Orgn] AND APY[Gene]") 
record = Entrez.read(handle) 
#ID
ID=record['IdList']
print(ID)
handle2 = Entrez.efetch(db="protein", id=ID[0], rettype='gb',retmode="txt")
record2 = SeqIO.read(handle2, "genbank")
handle.close()
print(record2.annotations['accessions'])
# Guarda info para otras db:
db_source = record2.annotations['db_source']
apy_prot = record2.annotations['accessions']


#ExPASy: herramienta para buscar y obtener:
from Bio import ExPASy 
from Bio import SwissProt

handle = ExPASy.get_sprot_raw(apy_prot[0])
record = SwissProt.read(handle)
print(record.__dict__.keys())
print (record.sequence_length)

#Entry name/locus También podemos acceder con el entry name o el locus

# Información extraída con efetch 
db_source.split(';')[0]
handle = ExPASy.get_sprot_raw('APY_AEDAE')
record = SwissProt.read(handle)
print (record.entry_name)
print (record.sequence_length)
# Si ha sido revisado
print (record.data_class)
print (record.organism)
print (record.sequence[:10]) # String, NO Bio.Seq



#Ejercicio 2:
#ExpASy obtiene un archivo crudo sin parsear
# SwissProt crea el objeto e¡con el archivo crudo
handle = ExPASy.get_sprot_raw('APY_AEDAE')
record = SwissProt.read(handle)
print(record.created)
print(record.annotation_update)
print(record.taxonomy_id)




#Con un file descargado en handle es un open al file.txt con 1 solo record
# Si el file tiene multiples records separados po // el handle es open
# y el archivo para leer un parse

#Crear objeto SeqRecord
import Bio.SeqRecord, Bio.Seq
# Se deben dar el id y seq
seqRec = Bio.SeqRecord.SeqRecord(seq=Bio.Seq.Seq(record.sequence),
                                 id=record.entry_name,
                                 name=record.organism,
                                 description=record.description)
print(seqRec.__dict__)


#Ejercicio 3:
seqRec = Bio.SeqRecord.SeqRecord(seq=Bio.Seq.Seq(record.sequence),
                                 id=record.entry_name,
                                 name=record.organism,
                                 description=record.description)

print(seqRec.format('fasta'))



# Prosite: db de familias, dominios
from Bio.ExPASy import Prosite
#Info de prosite en el file Swissprot
for reference in record.cross_references:
    if 'PROSITE' in reference:
        print(reference)

handle = ExPASy.get_prosite_raw("PS00785")
record = Prosite.read(handle)
print (record.name)

#Pedir documento para ver toda la info que da prosite
from Bio.ExPASy import Prodoc
handle = ExPASy.get_prosite_raw(documento)
record = Prodoc.read(handle)
print(record.text.replace('\\n','\n'))



# Scan prosite
sequence = "MEHKEVVLLLLLFLKSGQGEPLDDYVNTQGASLFSVTKKQLGAGSIEECAAKCEEDEEFTCRAFQYHSKEQQCVIMAENRKSSIIIRMRDVVLFEKKVYLSECKTGNGKNYRGTMSKTKN"
from Bio.ExPASy import ScanProsite
handle = ScanProsite.scan(seq=sequence)
result = ScanProsite.read(handle)
type(result)
print(result[0])

record.text








