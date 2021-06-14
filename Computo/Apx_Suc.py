# Daianna Gonzalez Padilla.
# Computo cientifico, 2do Semestre 2021, LCG: Proyecto Final
# Programa 1: Metodo de aproximaciones sucesivas para aproximar raices de una funcion

# Diccionario de variables:

# empezar permite reiniciar el programa
# fx y gx guardan las funciones dadas como input, f y g toma tal input y lo ve como una expresion en terminos de x
# k es el numero de decimales exactos que se desea obtener en la aproximacion
# n es el limite inferior del intervalo de busqueda y m el superior
# j es el incremento para la tabulacion
# tab es un diccionario que tiene como keys los x's de la tabulacion y a cada f(x) como values
# rangos es una lista que contiene a aquellos intervalos donde existe al menos una raiz
# i sirve como contador
# x,y,z guardan aproximaciones sucesivas Xn-1, X, Xn+1
# converg es una bandera que permite determinar si el metodo diverge o converge con el despeje dado en cierto intervalo
# xi y y son variables locales de la funcion print_graphs para graficar la funcion tabulando en un intervalo dado


print("Programa para aproximar raices de una funcion por aproximaciones sucesivas \n")
print("Consideraciones:\n")
print("*Considere que dentro de un intervalo dado puede no haber raices o haber mas de una\n")
print("*Con el incremento para tabular pueden no reconocerse raices incluso cuando existan, considere otros\n")
print("*Introduzca la funcion considerando la jerarquia de operaciones y las expresiones validas tales como pi, "
      "exp(), sin, arcsin\n")
print("*Con algunos despejes puede que al evaluar algunos valores de X no sean reales, si bien el metodo converge\n")
print("*El criterio de convergencia se basa en la reduccion del error con cada aproximacion y considera que los valores "
      "esten en el rango\n")
print("*Existen valores de x fuera del dominio de la funcion. Asegurese de no considerarlos dentro del intervalo\n")
print("*El error relativo mostrado es para la aproximacion final dada por el programa\n")

import sympy as sp
from sympy import *
import numpy as np
import matplotlib.pyplot as plt

# Cada que se oprima 1 se empieza el programa
empezar = int(input("Para comenzar oprima 1\n"))
while empezar == 1:

    # Funcion que imprime las graficas de la funcion y la raiz obtenida
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

    # Se pide el despeje g(x)
    def despeje():
        gx = input("Inserte el despeje de la variable x, asegurese de escribir la expresion correcta: ")
        return (simplify(gx))


    # Funcion para obtener la raiz en cada intervalo; se evalua que el metodo converga en ese intervalo para ese despeje
    def soluciones(intervalo, g, k):
        # Se obtienen las primeras 3 aproximaciones para evaluar si el error se reduce o no, suponiendo inicialmente
        # que converge (converg=True)
        x = round((intervalo[0] + intervalo[1]) / 2, k + 1)
        y = round(g.evalf(subs={'x': x}), k + 1)
        z = round(g.evalf(subs={'x': y}), k + 1)
        converg = True
        # Condicion de finalizacion de las iteraciones: cuando los dos ultimos X's posean los mismos k decimales,
        # se imprime la raiz y la grafica  g(x)
        while not int(z * (10 ** k)) == int(y * (10 ** k)):
            # Si |X(n+1)-X(n)| > |X(n)-X(n-1)|, para alguna iteracion, entonces el error se amplifica y el metodo diverge
            # Se considera cuando la diferencia es igual pues si esto pasa Xo=X2=X4.. y X1=X3... entonces tampoco se
            # converge
            # Si se converge pero en otro intervalo (z esta fuera de el), se considera que diverge
            if abs(z - y) >= abs(y - x) or (z < intervalo[0] or z > intervalo[1]):
                converg = False
                break
            x = y
            y = z
            # Se calcula la siguiente aproximacion
            z = round(g.evalf(subs={'x': y}), k + 1)

        # Si el metodo converge se imprime la solucion, el error relativo de la aproximacion final
        # y se grafica la solucion
        if converg == True:
            print("Raiz en el intervalo ", round(intervalo[0], 3), ":", round(intervalo[1], 3), "= ", round(z, k))
            print("Error relativo: ", np.format_float_scientific(abs((z - y) / z), precision=3, exp_digits=1))
            graficas = print_graphs(f, z)
        # Si no se converge se pide otro despeje u otro intervalo/incremento
        else:
            print("No converge en el intervalo ", round(intervalo[0], 3), ":", round(intervalo[1], 3))
            print("Oprima 1 para ingresar otro despeje, o 2 para otro intervalo/incremento\n")
            opcion = int(input())
            if opcion == 1:
                g = despeje()
                raices = soluciones(intervalo, g, k)
            else:
                inter_2 = main(f, k)
                raices = soluciones(inter_2, g, k)


    # Funcion principal: en esta se hace la tabulacion para todo X dentro del intervalo dado y con el incremento j
    def main(f, k):
        tab = {}
        rangos = []
        # f=simplify(f)
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

        # Se encuentran los intervalos i, i+j para los cuales hay un cambio de signo en f(x)
        # y se guardan en ranges para buscar la raiz en ellos
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
                inter_2 = main(f, k)
        # Si se encontraron raices, se buscan con la funcion soluciones para cada intervalo, con el despeje dado
        else:
            g = despeje()
            for intervalo in rangos:
                soluciones(intervalo, g, k)


    # Se pide la funcion f(x), se muestra la curva de f(x) para que el usuario determine el intervalo de busqueda de
    # la o las raices
    fx = input("Inserte la funcion, asegurese de usar las expresiones correctas y permitidas: \n")
    x = sp.symbols('x')
    f = simplify(fx)
    # Se pide el numero de decimales exactos
    k = int(input("Numero de decimales exactos: "))

    # Se imprime la grafica
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
    # Si la funcion no es continua en -10 a 10 no se imprime
    except TypeError:
        print("Considere que la funcion posee intervalos no continuos dentro de -10 y 10")

    # Se llama a la funcion principal para ingresar el intervalo y obtener las soluciones
    busqueda = main(f, k)

    empezar = int(input("Para ingresar otra funcion oprima 1, para salir 2\n"))
print("Gracias por usar este programa")