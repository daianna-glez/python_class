file=open("hello.py")
file_contents=file.read()

file.close()

# Agregar el nombre del archivo como variable
my_file_name = "../docs/dna.txt"
# Leer el archivo
my_file = open(my_file_name)
my_dna = my_file.read()
# Calular la longitud
dna_length = len(my_dna)
# Mostrar mensaje
print("The sequence is " + my_dna +  " and length is " + str(dna_length))
my_file.close()