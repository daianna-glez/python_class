# Daianna Gonzalez Padilla.
# Computo cientifico, 2do Semestre 2021, LCG: Proyecto Final
# Programa 2: Metodo de Newton-Raphson para aproximar raices de una funcion

# Diccionario de variables:

# empezar permite reiniciar el programa
# opcion permite elegir entre ingresar otra funcion o salir segun el valor que tome
# fx guarda las funcion dada como input, f toma tal input y lo ve como una expresion en terminos de x
# expr es la formula para la n+1-esima expresion
# k es el numero de decimales exactos que se desea obtener en la aproximacion
# n es el limite inferior del intervalo de busqueda y m el superior
# j es el incremento para la tabulacion
# tab es un diccionario que tiene como keys los x's de la tabulacion y a cada f(x) como values
# rangos es una lista que contiene a aquellos intervalos donde existe al menos una raiz
# i sirve como contador
# x,y son aproximaciones sucesivas  X, Xn+1
# xi y y son variables locales de la funcion print_graphs para graficar la funcion tabulando en un intervalo dado

import sympy as sp
from sympy import *
import numpy as np
import matplotlib.pyplot as plt

# Cada que se oprima 1 se empieza el programa
empezar = int(input("Para comenzar oprima 1\n"))
while empezar == 1:

    # Funcion que imprime la grafica de la funcion y la raiz obtenida
    def print_graphs(f, z):
        try:
            xi = np.linspace(-20, 20, 1000)
            y = [f.evalf(subs={'x': i}) for i in xi]
            plt.plot(xi, y)
            plt.axhline(0, color="black")
            plt.axvline(0, color="black")
            # Se definen los limites de la grafica mostrada
            plt.xlim(int(z) - 8, int(z) + 8)
            plt.ylim(-8, 8)
            # Se muestra el punto de la raiz
            plt.scatter([z, ], [f.evalf(subs={'x': z}), ], 50, color='red')
            plt.annotate(r'Raiz aproximada', xy=(z, 0), fontsize=8,
                         arrowprops=dict(arrowstyle="->"))
            print('Grafica de las funcion con la solucion obtenida:\n')
            plt.show()
        # Si la funcion no es continua en el intervalo cercano a la raiz no se imprime la grafica
        except TypeError:
            pass




    # Funcion para obtener la raiz en cada intervalo: se obtiene Xo y los siguientes X evaluando el valor previo de X
    # en expr
    def soluciones(intervalo, expr, f, k):
        x = round((intervalo[0] + intervalo[1]) / 2, k + 1)
        y = round(expr.evalf(subs={'x': x}), k + 1)
        # Condicion de finalizacion de las iteraciones: cuando los dos ultimos x posean los mismos k decimales
        while not int(y * (10 ** k)) == int(x * (10 ** k)):
            x = y
            y = round(expr.evalf(subs={'x': x}), k + 1)
        # Se imprime la raiz, el error relativo y la grafica con la solucion
        print("Raiz en el intervalo ", round(intervalo[0], 3), ":", round(intervalo[1], 3), "= ", round(y, k))
        print("Error relativo: ",np.format_float_scientific(abs((y - x) / y),precision = 3, exp_digits=3))
        graficas = print_graphs(f, y)


    # Funcion principal: en esta se hace la tabulacion para todo X dentro del intervalo dado y con el incremento j
    def main(f, expr, k):
        tab = {}
        rangos = []
        expr = simplify(expr)
        # Intervalo para tabular
        print("Intervalo:\n")
        n = float(input("min: "))
        m = float(input("max: "))

        # Tabulacion de valores para los x en el intervalo con incremento de j
        j = float(input("Ingrese el incremento para la tabulacion: "))
        i = n
        while i < m + j:
            tab[i] = round(f.evalf(subs={'x': i}), k + 1)
            i = i + j


        # Se encuentran los intervalos i, i+j para los cuales hay un cambio de signo en f(x) y se guardan en ranges
        # para buscar la raiz en ellos
        i = n
        while i < m:
            if ((tab[i] * tab[i + j]) < 0) and (i, i + j) not in rangos:
                rangos.append((i, i + j))
            i = i + j

        # Si no se encontraron cambios de signo en ningun intervalo, ranges esta vacio y no se encontraran raices,
        # se da otro intervalo u otra funcion
        if len(rangos) == 0:
            print(
                "Sin raices encontradas, oprima 1 para intentar con otro intervalo o incremento, o 2 para otra funcion/salir")
            opcion = int(input())
            if opcion == 1:
                inter_2 = main(f, expr, k)
        # Si se encontraron raices, se buscan con la funcion soluciones para cada intervalo, con el despeje dado
        else:
            for intervalo in rangos:
                soluciones(intervalo, expr,f, k)


    # Se pide la funcion f(x), se muestra la curva de f(x) para que el usuario determine el intervalo de busqueda de la o las raices
    fx = input("Inserte la funcion, asegurese de usar las expresiones correctas y permitidas: \n")
    x = sp.symbols('x')
    f = simplify(fx)
    # Se pide el numero de decimales exactos
    k = int(input("Numero de decimales exactos: "))
    #Se muestra la grafica de la funcion
    try:
        xi = np.linspace(-10, 10, 1000)
        y = [f.evalf(subs={'x': i}) for i in xi]
        plt.plot(xi, y)
        # Se definen los limites de la grafica mostrada
        plt.axhline(0, color="black")
        plt.axvline(0, color="black")
        plt.xlim(-10, 10)
        plt.ylim(-30, 30)
        print(
            "Con base en la grafica de la funcion determine el intervalo que contiene a la o las raices de interes.\n")
        plt.show()
    # Si la funcion no es continua en el intervalo de -10 a 10 no se imprime la grafica
    except TypeError:
        print("Considere que la funcion posee intervalos no continuos dentro de -10 y 10")

    # El valor del siguiente x se obtiene al evaluar el x anterior en expr
    expr = x - (f / f.diff())
    #Se llama a la funcion principal y se piden el intervalo de busqueda de la raiz
    busqueda = main(f, expr, k)

    empezar = int(input("Para ingresar otra funcion oprima 1, para salir 2\n"))
print("Gracias por usar este programa")