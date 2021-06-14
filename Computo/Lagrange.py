# Daianna Gonzalez Padilla.
# Computo cientifico, 2do Semestre 2021, LCG: Proyecto Final
# Programa 4: Metodo de interpolacion de Lagrange

# Diccionario de variables:

# empezar permite reiniciar el programa
# f guarda la funcion del polinomio y p a la expresion en terminos de x
# pos es un diccioanrio que contiene las coordenadas dadas por el usuario
# n es el numero de puntos dados por el usuario
# axis_x es una lista con las coordenadas en x de los puntos
# axis_y es una lista con las coordenadas en y de los puntos
# p1 es el producto de los (X-Xj) y p2 de (Xi-Xj) para i,j en {1,2,3,..,n} con i distinto de j
# i,j son contadores y representan a los indices de los X
# p es el polinimio obtenido de la interpolacion
# xi y y son variables para graficar la funcion tabulando en un intervalo dado
# X y Y son arreglos para graficar los puntos dados

print("Programa para interpolar una funcion por Lagrange \n")


import sympy as sp
from sympy import *
import pylab as pl
from pylab import *
import numpy as np


# Cada que se oprima 1 se empieza el programa
empezar = int(input("Para comenzar oprima 1\n"))
while empezar == 1:

    pos = {}
    axis_x = []
    axis_y = []
    n = int(input("A partir de cuantos puntos desea interpolar la curva?"))
    i = 0
    # Se piden y guardan las coordenadas de los puntos en las listas axis_x y axis_y, se guarda cada coordenada en pos
    while i < n:
        x = float(input("Ingrese la coordenada x para el punto " + str(i + 1) + " "))
        y = float(input("Ingrese la coordenada y para el punto " + str(i + 1) + " "))
        axis_x.append(x)
        axis_y.append(y)
        pos[i] = [x, y]
        i = i + 1


    # Funcion para evaluar un punto en el polinomio resultante
    def eval_punto(p):
        x = float(input("Ingrese el valor de x para el que quiere evaluar la funcion:\n"))
        print("El valor es: ", round(p.evalf(subs={'x': x}), 4))
        opcion = int(input("Para evaluar la funcion con otro x dado, oprima 1, si no 2:\n"))
        if opcion == 1:
            value = eval_punto(p)


    # f es el polinimio que describe a la curva que pasa por los puntos, es la suma de p1/p2*(Yi) con i en {1,2,3,...,n},

    # y f comienza en cero
    f = sp.Function('f')(x)
    f = 0
    # Se toma a cada coordenada Xi,Yi
    for i in pos.keys():
        # Los productos comienzan en 1: p1 es el producto de los (X-Xj) y p2 de (Xi-Xj) para i,j en {1,2,3,..,n} con
        # i distinto de j
        p1 = 1
        p2 = 1
        # Para cada coordenada se calula (X-Xj)/(Xi-Xj)*Yi, o sea p1/p2*Yi
        for j in pos.keys():
            # Se descarta cuando i=j pues de ser asi el denominador de la division p1/p2 es (Xi-Xj)=0 y la formula misma
            # no lo considera
            if not i == j:
                x = sp.symbols('x')
                # Se multiplica p1 por cada (X-Xj) y a p2 por cada (Xi-Xj)
                p1 = p1 * (x - axis_x[j])
                p2 = p2 * (axis_x[i] - axis_x[j])
        # Se suma a f cada p1/p2*Yi
        f = f + (p1 / p2 * (axis_y[i]))
    p = simplify(f)
    # Se imprime el polinomio resultante
    print("Polinomio: ", p)

    # Se muestra la grafica con los puntos dados y la curva resultante
    # Se acomoda el rango para que se vean todos los puntos y se evaluan en el polinomio obtenido p
    xi = np.linspace(int(min(axis_x)) - 2, int(max(axis_x)) + 2, 1000)
    y = [p.evalf(subs={'x': i}) for i in xi]
    plt.plot(xi, y)
    #Se imprimen los puntos originales dados
    X = np.array(axis_x)
    Y = np.array(axis_y)
    pl.scatter(X, Y, 50, color='red')
    pl.show()

    # Se pregunta si se quiere evaluar algun x en el polinomio dado
    opcion = int(
        input("Para evaluar la funcion con un  x dado, oprima 1, para interpolar con otros puntos/salir oprima 2:\n"))
    if opcion == 1:
        value = eval_punto(p)

    empezar = int(input("Para interpolar con otro conjunto de puntos oprima 1, para salir 2\n"))
print("Gracias por usar este programa")