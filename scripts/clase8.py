# Clase 8
## Usar matrices

import numpy as np

count_matrix = np.array([ [3, 3, 0],

                          [0, 0, 1],

                          [1, 1, 0],

                          [0, 0, 1],

                          [1, 0, 4],

                          [1, 2, 0]])
# Matriz de coexpresión: saber qué genes se expresan en cada celula (binaria)
# Donde la matriz sea mayor a 0, será 1, si no 0
count_matrix=np.where(count_matrix>0,1,0)
print(count_matrix)

# Comparar gen vs gen: producto punto
# Valores del gen 1 (hacia abajo)
print(count_matrix.T[0])
print(np.dot(count_matrix.T[0], count_matrix.T[0])
# Apilar:
print(np.vstack(count_matrix.T[0],count_matrix[0]))

# Suma de intersección de vectores:
print(sum((count_matrix[0]) & (count_matrix[0]))))


# producto punto de los dos vecotores
print("producto punto gen1vsgen1:", np.dot(count_matrix.T[0], count_matrix.T[0]))


#### Gen 1 vs gen 3
print(np.vstack((count_matrix.T[0],count_matrix.T[2])) )
# para visualizar
print("producto punto gen1 vs gen3:", np.dot(count_matrix.T[0], count_matrix.T[2]))
expresion=np.matmul(count_matrix.T, count_matrix)

import matplotlib.pyplot as plt
#Plot de matriz de coexpresión:
plt.inshow(expresion)
plt.colorbar()
plt.show()

import networkx as nx

G = nx.DiGraph(expresion)

nx.draw(G, node_size = 900,  with_labels=True)
plt.show()

#Seeds:
# Random para el seed de numpy
# Cada que se haga algo random se llama al seed
np.random.seed(10)
print(np.random.rand(1),np.random.rand(1), np.random.rand(1))



