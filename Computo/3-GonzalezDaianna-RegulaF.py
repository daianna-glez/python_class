# Daianna Gonzalez Padilla.
# Computo cientifico, 2do Semestre 2021, LCG: Proyecto Final
# Programa 3: Metodo de Regula Falsi para aproximar raices de una funcion

# Diccionario de variables:

# empezar permite reiniciar el programa
# fx guarda la funcion dada como input, f toma tal input y lo ve como una expresion en terminos de x
# k es el numero de decimales exactos que se desea obtener en la aproximacion
# n es el limite inferior del intervalo de busqueda y m el superior
# j es el incremento para la tabulacion
# tab es un diccionario que tiene como keys los x's de la tabulacion y a cada f(x) como values
# rangos es una lista que contiene a aquellos intervalos donde existe al menos una raiz
# i sirve como contador
# x,y son aproximaciones sucesivas Xn-1, X, Xn+1
# q es el pivote con el cual se obtendra la solucion
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


    # Funcion para obtener la raiz en cada intervalo: se obtiene Xo y los siguientes X
    def soluciones(intervalo, f, k):
        # Se determina el valor de Xo dado por a*f(b)-b*f(a)/f(b)-f(a)
        a = intervalo[0]
        b = intervalo[1]
        x = ((a * f.evalf(subs={'x': b})) - (b * f.evalf(subs={'x': a}))) / (
                    (f.evalf(subs={'x': b})) - (f.evalf(subs={'x': a})))
        # Se determina que extremo del intervalo sera el pivote
        if f.evalf(subs={'x': intervalo[0]}) * f.evalf(subs={'x': x}) < 0:
            q = intervalo[0]
        elif f.evalf(subs={'x': intervalo[1]}) * f.evalf(subs={'x': x}) < 0:
            q = intervalo[1]
        # El siguiente X (denotado por y) esta dado por x*f(k)-k*f(x)/f(k)-f(x)
        y = round((((x * f.evalf(subs={'x': q})) - (q * f.evalf(subs={'x': x}))) / (
                    (f.evalf(subs={'x': q})) - (f.evalf(subs={'x': x})))), k + 1)
        # Condicion de finalizacion de las iteraciones: cuando los dos ultimos x posean los mismos k decimales
        while not int(y * (10 ** k)) == int(x * (10 ** k)):
            # Si F(Xn)*F(Xn+1)<0, se cambia el pivote y ahora es Xn
            if f.evalf(subs={'x': x}) * f.evalf(subs={'x': y}) < 0:
                q = x
            x = y
            y = round((((x * f.evalf(subs={'x': q})) - (q * f.evalf(subs={'x': x}))) / (
                        (f.evalf(subs={'x': q})) - (f.evalf(subs={'x': x})))), k + 1)
        #Se imprime la raiz, el error relativo y la grafica con la solucion
        print("Raiz en el intervalo ", round(intervalo[0], 3), ":", round(intervalo[1], 3), "= ", round(y, k))
        print("Error relativo: ", np.format_float_scientific(abs((y - x) / y), precision=1, exp_digits=2))
        graficas = print_graphs(f, y)


    # Funcion principal: en esta se hace la tabulacion para todo X dentro del intervalo dado y con el incremento dado
    def main(f, k):
        tab = {}
        rangos = []
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

        # Se encuentran los intervalos i, i+j para los cuales hay un cambio de signo en f(x) y se guardan en ranges para buscar la raiz en ellos
        i = n
        while i < m:
            if ((tab[i] * tab[i + j]) < 0) and (i, i + j) not in rangos:
                rangos.append((i, i + j))
            i = i + j

        # Si no se encontraron cambios de signo en ningun intervalo, ranges esta vacio y no se encontraran raices, se da otro intervalo u otra funcion
        if len(rangos) == 0:
            print(
                "Sin raices encontradas, oprima 1 para intentar con otro intervalo o incremento, o 2 para otra funcion/salir")
            opcion = int(input())
            if opcion == 1:
                inter_2 = main(f, k)
        # Si se encontraron raices, se buscan con la funcion soluciones para cada intervalo, con el despeje dado
        else:
            for intervalo in rangos:
                soluciones(intervalo, f, k)


    # Se pide la funcion f(x), se muestra la curva de f(x) para que el usuario determine el intervalo de busqueda de la o las raices
    fx = input("Inserte la funcion, asegurese de usar las expresiones correctas y permitidas: \n")
    x = sp.symbols('x')
    f = simplify(fx)
    # Se pide el numero de decimales exactos
    k = int(input("Numero de decimales exactos: "))
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


    # Se llama a la funcion busqueda
    busqueda = main(f, k)

    empezar = int(input("Para ingresar otra funcion oprima 1, para salir 2\n"))
print("Gracias por usar este programa")