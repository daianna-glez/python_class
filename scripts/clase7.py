# Clase 7: NumPy
# NumPy maneja arrays
## Shape: da el tamaño del arreglo hay 3 ejes: el 0,1,2... segun filas, columnas, 
'''
shape: (x,y,z)
No se mecesita iterar, son objetos 
'''
import numpy as np
# Al declarar array debe estar en un solo objeto con listas o tuplas que encierren
# a los elementos del array
array_1D=np.array([1,2,3])
array_2D=np.array(([1,2,3],[4,5,6])) # Guarda una matriz
print(array_1D,array_2D)

# Preguntar no de dimensiones: ndim
# Longitud de las filas: len()
# Forma del arreglo: shape

# biomasa en unidades de absorbancia
ecoli_matraz = np.array([0.1, 0.15, 0.19, 0.5,  0.9, 1.4, 1.8, 2.1, 2.3])
print(ecoli_matraz.ndim, ecoli_matraz.shape, len(ecoli_matraz))


# Arreglo 2D
array_2D=np.array([(1,2,3),(4,5,6)])
print(array_2D.ndim, array_2D.shape, len(array_2D))

# Arreglo 3D
# Biomasa en unidades de absorbancia  (OD600)
ecoli_m_b = np.array([[0.1, 0.15, 0.19, 0.5,  # Matraz 250 mL
                       0.9, 1.4, 1.8, 2.1, 2.3],
                      [0.1, 0.17, 0.2, 0.53,  # Biorreactor 50 L
                       0.97, 1.43, 1.8, 2.1,  2.8],
                      [0.1, 0.17, 0.2, 0.52,  # B. alimentado 50 L
                       0.95, 1.41, 1.8, 2.2,  2.8]
                    ])
print(ecoli_m_b.ndim, ecoli_m_b.shape,len(ecoli_m_b))

# Operaciones con un único valor y entre matrices:
#Multiplicar cada elemeto por un valor
ecoli_matraz_gL = ecoli_matraz*0.39
print(ecoli_matraz_gL)

# Ejercicio 1:
produccion = np.array([[16, 14], [12, 9]])
total = produccion*2
print(produccion+produccion)  # Por bacteria y metabolito
# np.sum para sumar la matriz con respecto a un eje
print(np.sum(produccion*2, axis= 0)) # Por metabolito
print(np.sum(produccion*2, axis= 1)) # Por bacteria
contaminado = produccion/2
total_real = total - contaminado


# Ejercicio 2:
consumo =  np.array([[7, 3],[5, 2]])
sobra=total_real-consumo
consumo =  np.array([[7, 3],[7, 3]])
total_real - consumo # Por bacteria

# Reciclaje: broadcasting 
consumo =  np.array([7, 3])
sobra=total_real - consumo

# Otras operaciones:
## Potencia
print(total_real**2)
## Transpuesta del array
print(total_real.T)
## Suma de todo:
print(total_real.sum())
## Min y max:
print( total_real.max(), np.max(total_real) )
## Exp y raiz cuadrada
print(np.exp(total_real))
print(np.sqrt(total_real))
## Trigonométricas
print(np.sin(np.array([np.pi, np.pi/2])))
print(np.arcsin(np.array([0.0, 1.0])))
## Calcular y asignar
## Redondeo
redondear= np.array([1.1, 1.5, 1.9, 2.5])
print(np.floor(redondear))
print(np.ceil(redondear))
print(np.round(redondear))

#EJERCICIO 1:
produccion = np.array([ [5,3], [11,7], [4,9], [2,6]])
costos = np.array([3.5, 5, 7, 4.3])
sobra=costos/produccion.T
print(sobra.T)
print(np.min(sobra.T))

# Tipos de datos:
from sys import getsizeof
# Decimales:
np_float = np.array([1.0, 2.0, 3.0, 4.0])
print("Tipo de dato\t", np_float.dtype, 
"\nTamaño en bytes\t", getsizeof(np_float))
# Enteros:
np_int = np.array([1, 2, 3, 4])
print("Tipo de dato\t", np_int.dtype, 
"\nTamaño en bytes\t", getsizeof(np_int))
# También podemos especificar el tipo de dato
np_float = np.array([1, 2, 3, 4],dtype='float64')
print("Tipo de dato\t", np_float.dtype, 
"\nTamaño en bytes\t", getsizeof(np_float))

#Tipo de dato booleano:
bool_np = np.array([True, False, True, False]) 
bool_np.dtype
## Para acceder a los Trues:
np_int=np.array([1, 2, 3, 4])
print(np_int[bool_np])
# Obtener array de Trues y Falses <- vector booleano
new=np_int <3
print(new)
# Obtener los elementos True
print(np_int[np_int <3])
print(np_int[np_int == np_int.max()])
print(np_int[(np_int <2 ) | (np_int >3)])


# EJERICIO 2:
sobra=sobra.T
print(sobra[(sobra==sobra.min())| (sobra==sobra.max())])

# Tipo de dato complejo
num_1 = np.array([3+6j])   
num_2 = np.array([7+2j])   
print(num_1.dtype)
# Parte real:
print(num_1.real)
# Parte imaginaria:
print(num_1.imag)

# Tipo de dato fecha:
dias = np.datetime64('2005-02-25')

# Acceder al array ID: como lista o matriz
# lista:
print(ecoli_matraz[2])
print(ecoli_matraz[2:5])
print(ecoli_matraz[0:6:2]) # Del 0 al 6 de 2 en 2
print(ecoli_matraz[:6:2]) # Del 0 al 6 de 2 en 2
print(ecoli_matraz[::-1]) # De fin a inicio
print(ecoli_matraz[::2]) # de inicio a fin de 2 en 2

# Acceder al arreglo 2D:
print(produccion[2],'\n', produccion[2:4])
print(produccion[0:6:2])  # Del 0 al 6 de 2 en 2
# Indices por coordenadas columna, fila:
print(produccion[2][1]) # o con produccion[2,1]
#produccion[0:6:2][1] (accede a las filas, luego se queda con la segunda)
# y no es produccion[0:6:2, 1] que accede a las filas y se queda con sus segundos
# elementos

#
print(produccion[-2:,-2:])
# Saca las filas de 2 en 2 y luego accede a las columnas 1 y sigue de 2-2
print(produccion[::2,1::2])


#... <- para cosas muy grandes como 10D
a_3D = np.array([[[  1,  2,  3],
               [ 11, 12, 13]],
              [[101, 102, 103],
               [1001, 1002, 1003]]])
print(a_3D[1, ...])  # a_3D[1, :, :] o a_3D[1] <- 2da rebanada
print(a_3D[0,...])
print(a_3D[ ...,2]  # Elementos de las ultimas columnas

# Crear arrays
a=np.arange(10)
a=np.arange(0, 10, 2) <- sin tocar al 10 con paso de 2-2
# Inicio, fin y en cuantos elementos partirlo con linespace
b=np.linspace(0, 8, 5)
#Problemas con redondeos <-np.arange(0.5, 0.8, 0.1) llega al 0.8
# => usar un numero al que no llegue como 0.75 o usar linespace

# Random:
print(np.random.randint(0, 10, 3)) # sin incluir al 10
# random.normal: no.valores, media, ds

# Arrays vs listas
## Arrays son homogeneos, lista y df heterogeneas
## No se puede aceder por coords en una lista como [1,2]
## En listas []*2 <- las duplica (son copias) y se modifican juntas


## Array estructurado: generar objetos struct
mascotas = np.array([('Freya', 6, 6.5), ('Senna', 1, 2.5)],
       dtype=[('nombre', (np.str_, 10)), ('edad', np.int32), ('peso', np.float64)])
# Establecer tipo de dato con dtype para los 1eros elementos, 2dos, etc
## Se indica nombre del tipo de dato, tamaño
# Ahora se puede ordenar el array
sort_age = np.sort(mascotas, order='edad')
sort_name = np.sort(mascotas, order='nombre')


# Tarea: arrays estructurados de cada tabla y los costos individuales por gen y temp
      







