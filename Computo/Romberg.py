# Daianna Gonzalez Padilla.
# Computo cientifico, 2do Semestre 2021, LCG: Proyecto Final
# Programa 8: Metodo de integracion por el metodo de Romberg

# Diccionario de variables:

# empezar permite reiniciar el programa
# fx guarda las funcion dada como input, f toma tal input y lo ve como una expresion en terminos de x
# a es el limite inferior del intervalo de busqueda y b el superior
# d es el numero de decimales de exactitud
# lista es una lista que contiene a las aproximaciones que se van obteniendo
# c es la constante que multiplica a cada jk para obtener las aproximaciones Ik
# los jk son los valores por los cuales se multiplica c para obtener los Ik
# k, contador, p,q son contadores
# s es un acumulador
# m es el log2(k) con k potencia de 2
# i1 son los Ik
# i2 son las refinaciones de los Ik
# xi y y permiten graficar la funcion tabulando en un intervalo dado
# X y Y son arreglos para mostrar el area bajo la curva en el intervalo a-b



import sympy as sp
from sympy import *
import math
import numpy as np
import matplotlib.pyplot as plt


print("Programa para aproximar la integral definida de una funcion por Romberg\n")

# Cada que se oprima 1 se empieza el programa
empezar = int(input("Para comenzar oprima 1\n"))
while empezar == 1:

    # Se pide la funcion f(x) que se quiere integrar y el intervalo a-b para ello, asi como el numero de decimales exactos
    fx = input("Inserte la funcion que quiere integrar, asegurese de usar las expresiones correctas y permitidas: \n")
    x = sp.symbols('x')
    f = simplify(fx)
    a = float(input("Ingrese el limite inferior del intervalo:\n"))
    b = float(input("Ingrese el limite superior  del intervalo:\n"))
    d = int(input("Ingrese el numero de decimales exactos:\n"))

    # En lista se iran guardando las aproximaciones obtenidas
    lista = []
    # Se calcula c=b-a
    c = b - a

    # Se calcula el primer Ik dado por Ik=c(jk) con jk=1/2*(f(a) + f(b))
    j = (1 / 2) * (f.evalf(subs={'x': a}) + f.evalf(subs={'x': b}))
    i1 = c * j
    lista.append(i1)

    # Se determinan los Ik siguientes desde la segunda iteracion
    k = 2
    contador = 2
    # Condicion de finalizacion: cuando In e In+1 tengan los mismos d decimales, esto se cumple cuando paro es verdadero
    paro = False
    while paro == False:
        s = 0
        p = 1
        # Se calcula j para la iteracion k; para k>1 jk= j(k-1)+ suma de f(a + c*m/k), donde m son los impares del 1 al k-1
        while p < k:
            s = s + (f.evalf(subs={'x': a + (c * (p / k))}))
            p = p + 2
        j = j + s
        # Nuevamente Ik=(1/k)c(jk)
        i1 = (1 / k) * c * j
        lista.append(i1)

        # Se calculan las refinaciones de Ik, para cado uno de ellos hay log2(k) refinaciones (I*,I**, etc.)
        p = 2
        q = 1
        m = math.log(k, 2)
        while q <= m:
            # Cada refinacion I2 es (1/(2**p)-1)*(4I1 - Ic) con p par empezando en 2, I1 es la aproximacion inmediata
            # anterior e Ic es la aproximacion obtenida hace 2,3,4 ... aproximaciones, dada por la iteracion en la que se este
            i2 = (1 / ((2 ** p) - 1)) * (((2 ** p) * (i1)) - lista[-contador])
            lista.append(i2)
            #Se evalua la condicion de paro
            if int(i1 * (10 ** d)) == int(i2 * (10 ** d)):
                paro = True
                break
            i1 = i2
            p = p + 2
            q = q + 1

        # El contador incrementa con cada iteracion y k son potencias de 2
        contador = contador + 1
        k = k * 2

    # Se imprime la aproximacion de la integral en ese intervalo para la funcion dada
    print("Aproximacion: ", round(lista[-1], d))

    # Se imprime la grafica de la integral de la funcion en tal intervalo
    # Se evaluan valores de x en f para graficar la funcion
    try:
        xi = np.linspace(a - 1, b + 1, 1000)
        y = [f.evalf(subs={'x': i}) for i in xi]
        plt.plot(xi, y)
        plt.axhline(0, color="black")
        plt.axvline(0, color="black")
        # Arreglos con los puntos en la curva que delimitan el area bajo la curva
        X = np.arange(int(a), int(b) + 0.1, 0.01)
        X = np.array(X, dtype=float)
        Y = [f.evalf(subs={'x': i}) for i in X]
        Y = np.array(Y, dtype=float)
        # Se definen los limites de la grafica mostrada
        plt.xlim(int(a) - 1, int(b) + 1)
        plt.ylim(float(min(Y)) - 1, float(max(Y)) + 1)
        # Se rellena el area bajo la curva en el intervalo de a -> b
        plt.fill_between(X, Y, color='r')
        plt.title('Integral definida de la funcion', fontsize=10)
        plt.show()
    # Si la funcion no es continua en el a-1, b+1, no se imprime la grafica
    except TypeError:
        pass

    empezar = int(input("Para ingresar otra funcion oprima 1, para salir 2\n"))
print("Gracias por usar este programa")