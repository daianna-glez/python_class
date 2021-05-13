'''
NAME
    drosophila_genespy

VERSION
    [1.0]

AUTHOR
    Daianna Gonzalez Padilla <daianna@lcg.unam.mx>

DESCRIPTION
    This programs takes a set of genes and gives a list of genes that fulfill certain requirements

CATEGORY
     Genes analysis

USAGE
    None

ARGUMENTS
    This program doesn't take arguments

INPUT
    The input text file, with each data separated by tab and whose columns are the name of the organism from
    which the gene is, the gene sequence, the gene name and its expression value.

OUTPUT
    List of genes with the specified characteristics

EXAMPLES
    Example 1: gets the following data:  Drosophila yakuba	cgcgcgctcgcgcatacggcctaatgcgcgcgctagcgatgc	hdt739	85
    and evaluate its AT content, the gene name, the organism name and the expression data.

SOURCE
    https://github.com/daianna21/python_class/blob/master/scripts/drosophila_genes.py


'''

#The file name is asked and all of its lines are saved
input_file_name = input("Insert the file name:\n")
input_file=open(input_file_name,"r")
all_lines = input_file.readlines()
input_file.close()


#Function to calculate AT content of a gene, takes the dna sequence and the figures to round
def get_at_content(dna, sig_figs):
    length = len(dna)
    a_count = dna.upper().count('A')
    t_count = dna.upper().count('T')
    at_content = (a_count + t_count) / length *100
    return round(at_content, sig_figs)


print("Genes that don not belong to Drosophila melanogaster or Drosophila simulans:\n")
for line in all_lines:
  #Create a list with the elements of each line
  line=line.split('\t')
  #Element 0 is the organism name
  if line[0]=="Drosophila melanogaster" or line[0]== "Drosophila simulans":
    #Print the gene
    print(line[2])

print("\nGenes of 90-110 nt:\n")
for line in all_lines:
  line=line.split('\t')
  #Element 1 is the dna sequence
  if len(line[1])>=90 and len(line[1])<=110:
    print(line[2])

print("\nGenes whose AT content is lower than 50 and whose expression level is higher than 200:\n")
for line in all_lines:
  line=line.split('\t')
  #Element 3 is the expression value of the gene
  if get_at_content(line[1],2)<50 and int(line[3])>200:
    print(line[2])

print("\nGenes whose name starts with k or h, except those belonging to Drosophila melanogaster:\n")
for line in all_lines:
  line=line.split('\t')
  #Evaluate if the gene name starts with those letters and is not part of a certain organism
  gene_name=line[2]
  if (gene_name.startswith('h') or gene_name.startswith('k') ) and not line[0]=="Drosophila melanogaster":
    print(line[2])

print("\nCualitative AT content of each gen:\n")
for line in all_lines:
  line=line.split('\t')
  #Depending on the AT content, the message
  if get_at_content(line[1],2)>65:
    print(line[2],": High")
  elif get_at_content(line[1],2)<45:
    print(line[2],": Low")
  else:
    print(line[2],": Medium")
