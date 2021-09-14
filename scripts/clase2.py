#Clase 2:
import Bio
from Bio import SeqIO
#Ejercicio 3: Leer archivos Fastq

mala_calidad=[]
umbral=40
promedio=0
for record in SeqIO.parse("C:/Users/hp/Downloads/sample.fastq", "fastq"):
    promedio=sum( record.letter_annotations["phred_quality"])/len(record.letter_annotations["phred_quality"])
    if (promedio<umbral):
        temp=(promedio,record.id)
        mala_calidad.append(temp)
print ('{} secuencias con promedio menor a umbral:{}\n'.format(len(mala_calidad),umbral))
