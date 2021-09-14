#Clase 1: Biopython, secuencias 

import Bio.Seq
from Bio.Seq import Seq
seqobj = Seq('ATGCGATCGAGCTAA')
seqobj
seq_str = str(seqobj)  # convertir con str()
print('{} tiene {} nucleotidos'.format( seq_str , len(seq_str)))

#Inmutabilidad:
print(seqobj[0])
# Al hacer seqobj[0]='T' marca error por la inmutabilidad
from Bio.Seq import MutableSeq

mutable=MutableSeq(seqobj)
mutable[0]='T'
print(mutable)

##Secuencia complementaria:
print(seqobj.complement())

##Secuencia reverso complementaria:
print(seqobj.reverse_complement())

##Secuencia traducida:
###(empieza a traducir desde inicio sin importar el codón de inicio)
####to_stop para en codón de paro, de lo contrario continua y escribe "*"
print(seqobj.translate(to_stop = True))
print(seqobj.translate(to_stop = False))

##Secuencia transcrita a RNA:
rna = seqobj.transcribe()
print(rna)

##Secuencia retrotranscrita
print(rna.back_transcribe())



#Extraer subsecuencias
##Extraer el primer codon
print(seqobj[0:3])

##Obtener los codones con expresiones regulares
import re
for codon in re.findall(r"(.{3})", str(seqobj)):
    print(codon)


    
#Buscar patrones:
    
from Bio import SeqUtils
## Patrón a buscar
pattern = Seq("ACG")
## Secuencia donde buscaremos patrón
sequence = Seq("ATGCGCGACGGCGTGATCAGCTTATAGCCGTACGACTGCTGC")
## Busqueda patrón en secuencia tipo STRING!!
results = SeqUtils.nt_search(str(sequence), pattern)
print(results)  # patrón y posiciones

##Buscar patron reverso complementario
# reverse_complement()
results_rc = SeqUtils.nt_search(str(sequence),
                pattern.reverse_complement())
print(results_rc)





#Contenido de GC
from Bio.SeqUtils import GC
print(GC(seqobj))
print(GC('AGCGTGCA'))



#calcular masa molecular
from Bio.SeqUtils import molecular_weight
print(molecular_weight(seqobj))

#Contar nt
##No sobrelapa, por lo que en "AAAA" dará 2 para "AA"
print(seqobj.count("ATG")) 
print(seqobj.find("GC"))  # nos da la primera posición que encuentra 
# (-1 si no encuentra)


#Ejercicio 1: Obtener cadena protéica de cualquiera de sus ORFs
secuencia=Seq('AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG')
secuencia_rev=secuencia.reverse_complement()
inicio=Seq('ATG')

#Regresa las posiciones de los codones de inicio
posicion=SeqUtils.nt_search(str(secuencia), inicio)
posicion_rev=SeqUtils.nt_search(str(secuencia_rev), inicio)

proteins=[]
#Empieza en 1 porque el elemento 0 de posisicones es el patron 
for i in range(1,len(posicion)):
    #Proteina desde cada codon de inicio hasta el de paro
    sec_prot=secuencia[i:]
    if len(sec_prot)%3==0:
        proteins.append(sec_prot.translate(to_stop=True))
#####REVISAR
for i in range(1,len(posicion_rev)):
    sec_prot=secuencia[i:]
    if len(sec_prot)%3==0:
        proteins.append(sec_prot.translate(to_stop=True))
print(proteins)



#Motifs
##Supongamos que tenemos las instancias de un DNA motif

from Bio import motifs
instances = [Seq("TACAA"),
            Seq("TACGC"),
            Seq("TACAC"),
            Seq("TACCC"),
            Seq("AACCC"),
            Seq("AATGC"),
            Seq("AATGC"),
           ]
##Crear un motif:
m = motifs.create(instances)
print(m)





#Objetos SeqRecord
from Bio import SeqIO
#dirección de archivo
filename = "C:/Users/hp/Downloads/seq.nt.fa"
#por cada record queremos ID, longitud seq y traducción
#Se pone el nombre del archivo y el formato en el que esta
for seq_record in SeqIO.parse(filename, "fasta"):
    print('ID {}'.format(seq_record.id))
    print('len {}'.format(len(seq_record)))
    print('Traducción {}'.format(seq_record.seq.translate(to_stop=False)))
#  (to_stop=False,cds=True)
#cds revisa que tenga un codón de inicio

###Secuencia es de tipo Seq por lo que podemos utilizar lo que vimos en Bio.Seq:
'''
seq_record.seq.complement()
seq_record.seq.reverse_complement()
rna = seq_record.seq.transcribe()
rna.back_transcribe()
'''


#Leer archivo en diccionario
#id_dict = SeqIO.to_dict(SeqIO.parse('archivos_trabajo/seq.nt.fa', 'fasta'))
#print(id_dict['seq3'])


#Ejercicio 2:
id_dict = SeqIO.to_dict(SeqIO.parse('C:/Users/hp/Downloads/seq.nt.fa', 'fasta'))
##Obtener los codones 
for i in id_dict:
    print('>',format(i))
    sec=id_dict[i].seq
    #Para los 3 marcos de lectura de la forward sec[1:] y sec[2:]
    for codon in re.findall(r"(.{3})", str(sec[0:len(sec)])):
        print(codon, end='\t')
    print('')



#FASTA con listas:
id_list = list(SeqIO.parse("C:/Users/hp/Downloads/seq.nt.fa", "fasta"))
print(id_list[-1].id, '\n')






#Ejercicio 3: Leer archivos Fastq

def corte_calidad(path, umbral):
    mala_calidad=[]
    for record in SeqIO.parse("C:/Users/hp/Downloads/sample.fastq", "fastq"):
       promedio=sum( record.letter_annotations["phred_quality"])/len(record.letter_annotations["phred_quality"])
       if (promedio<umbral):
           temp=(promedio,record.id)
           mala_calidad.append(temp)
    print ('{} secuencias con promedio menor a umbral:{}\n'.format(len(mala_calidad),umbral)


















